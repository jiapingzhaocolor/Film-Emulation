
#include "ofxsImageEffect.h"
#include "ofxsProcessing.h"
#include "ofxsMultiThread.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <array>

#include "luts_embedded.h"
#include "gamut_matrices.h"
#include "transfer_functions.h"

#define kPluginName "OpenDRT Film Pipeline"
#define kPluginGrouping "Color"
#define kPluginDescription "Input color-management into DaVinci Wide Gamut, then optional Negative / Color Separation / Print LUT stages."
#define kPluginIdentifier "com.opendrt.filmpipeline"
#define kPluginVersionMajor 1
#define kPluginVersionMinor 0

// Parameter names
static const char* kParamInGamut      = "in_gamut";
static const char* kParamInOETF       = "in_oetf";

static const char* kParamNegEnable    = "neg_enable";
static const char* kParamNegLut       = "neg_lut";
static const char* kParamNegBlend     = "neg_blend";

static const char* kParamSepEnable    = "sep_enable";
static const char* kParamSepStyle     = "sep_style";
static const char* kParamSepBlend     = "sep_blend";

static const char* kParamPrintEnable  = "print_enable";
static const char* kParamPrintLut     = "print_lut";
static const char* kParamPrintBlend   = "print_blend";

struct float3 { float x,y,z; };

static inline float3 make_float3(float x,float y,float z){ return {x,y,z}; }
static inline float3 operator+(const float3& a,const float3& b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline float3 operator-(const float3& a,const float3& b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline float3 operator*(const float3& a,float s){ return {a.x*s,a.y*s,a.z*s}; }
static inline float3 operator*(float s,const float3& a){ return a*s; }
static inline float3 operator/(const float3& a,float s){ return {a.x/s,a.y/s,a.z/s}; }

static inline float clampf(float v, float lo, float hi){ return std::min(std::max(v, lo), hi); }
static inline float3 clamp3(const float3& v, float lo, float hi){ return {clampf(v.x,lo,hi), clampf(v.y,lo,hi), clampf(v.z,lo,hi)}; }
static inline float3 lerp3(const float3& a, const float3& b, float t){ return a*(1.0f-t) + b*t; }
static inline float luma_rec709(const float3& v){ return 0.2126f*v.x + 0.7152f*v.y + 0.0722f*v.z; }

static inline float3 mat_apply(const Mat3& M, const float3& v){
    std::array<float,3> r = mat_vec(M, {v.x,v.y,v.z});
    return {r[0], r[1], r[2]};
}

// 3D LUT sampler (tetrahedral). LUT domain assumed [0,1].
struct Lut3D {
    const float* data; // flattened RGB triples, length = size^3*3
    int size;          // e.g., 33
};

static inline float3 lut_fetch(const Lut3D& lut, int r, int g, int b){
    const int N = lut.size;
    r = std::clamp(r, 0, N-1);
    g = std::clamp(g, 0, N-1);
    b = std::clamp(b, 0, N-1);
    // .cube order: blue fastest, then green, then red (common convention)
    const int idx = ((r * N + g) * N + b) * 3;
    return { lut.data[idx], lut.data[idx+1], lut.data[idx+2] };
}

static inline float3 lut_sample_tetra(const Lut3D& lut, const float3& in){
    const int N = lut.size;
    float3 x = clamp3(in, 0.0f, 1.0f);
    float fx = x.x * (N - 1);
    float fy = x.y * (N - 1);
    float fz = x.z * (N - 1);

    int ix = (int)std::floor(fx);
    int iy = (int)std::floor(fy);
    int iz = (int)std::floor(fz);

    float dx = fx - ix;
    float dy = fy - iy;
    float dz = fz - iz;

    // Corners
    float3 c000 = lut_fetch(lut, ix,   iy,   iz);
    float3 c100 = lut_fetch(lut, ix+1, iy,   iz);
    float3 c010 = lut_fetch(lut, ix,   iy+1, iz);
    float3 c001 = lut_fetch(lut, ix,   iy,   iz+1);
    float3 c110 = lut_fetch(lut, ix+1, iy+1, iz);
    float3 c101 = lut_fetch(lut, ix+1, iy,   iz+1);
    float3 c011 = lut_fetch(lut, ix,   iy+1, iz+1);
    float3 c111 = lut_fetch(lut, ix+1, iy+1, iz+1);

    // Tetrahedral interpolation (based on ordering of fractional parts)
    float3 out;
    if (dx >= dy) {
        if (dy >= dz) {
            // x >= y >= z
            out = c000
                + (c100 - c000) * dx
                + (c110 - c100) * dy
                + (c111 - c110) * dz;
        } else if (dx >= dz) {
            // x >= z > y
            out = c000
                + (c100 - c000) * dx
                + (c101 - c100) * dz
                + (c111 - c101) * dy;
        } else {
            // z > x >= y
            out = c000
                + (c001 - c000) * dz
                + (c101 - c001) * dx
                + (c111 - c101) * dy;
        }
    } else { // dy > dx
        if (dx >= dz) {
            // y > x >= z
            out = c000
                + (c010 - c000) * dy
                + (c110 - c010) * dx
                + (c111 - c110) * dz;
        } else if (dy >= dz) {
            // y >= z > x
            out = c000
                + (c010 - c000) * dy
                + (c011 - c010) * dz
                + (c111 - c011) * dx;
        } else {
            // z > y > x
            out = c000
                + (c001 - c000) * dz
                + (c011 - c001) * dy
                + (c111 - c011) * dx;
        }
    }
    return out;
}

// Convert input to DaVinciWG + Davinci Intermediate (encoded), using OpenDRT matrices + transfer decode.
static inline float3 input_to_dwg_intermediate(const float3& rgb_in, int inGamutIdx, int inOetfIdx){
    // Decode input transfer to linear
    float3 lin = {
        decode_input_oetf(inOetfIdx, rgb_in.x),
        decode_input_oetf(inOetfIdx, rgb_in.y),
        decode_input_oetf(inOetfIdx, rgb_in.z)
    };

    // Convert input gamut -> XYZ
    Mat3 M_in_to_xyz = inputGamutToXYZ(inGamutIdx);
    float3 xyz = mat_apply(M_in_to_xyz, lin);

    // Convert XYZ -> DaVinciWG (linear)
    Mat3 M_xyz_to_dwg = xyzToDaVinciWG();
    float3 dwg_lin = mat_apply(M_xyz_to_dwg, xyz);

    // Encode to Davinci Intermediate (0..~)
    float3 dwg_di = {
        encode_davinci_intermediate(dwg_lin.x),
        encode_davinci_intermediate(dwg_lin.y),
        encode_davinci_intermediate(dwg_lin.z)
    };
    return dwg_di;
}

// Baseline output: DaVinciWG+DI -> Rec709 gamma 2.4
static inline float3 dwg_di_to_rec709_24(const float3& dwg_di){
    // Decode DI to linear
    float3 dwg_lin = {
        decode_davinci_intermediate(dwg_di.x),
        decode_davinci_intermediate(dwg_di.y),
        decode_davinci_intermediate(dwg_di.z)
    };

    // DaVinciWG linear -> XYZ
    // We have DWG->XYZ matrix directly from OpenDRT (as input matrix)
    float3 xyz = mat_apply(matrix_davinciwg_to_xyz, dwg_lin);

    // XYZ -> Rec709 linear
    Mat3 M_xyz_to_rec709 = xyzToRec709();
    float3 rec_lin = mat_apply(M_xyz_to_rec709, xyz);

    // Encode Rec709 2.4
    float3 rec_24 = { encode_rec709_24(rec_lin.x), encode_rec709_24(rec_lin.y), encode_rec709_24(rec_lin.z) };
    return rec_24;
}

enum NegLutChoice { Neg_Cthulhu=0, Neg_Lilith=1, Neg_Tsathoggua=2, Neg_Yig=3 };
enum SepChoice { Sep_Hydra=0, Sep_Oorn=1, Sep_Zhar=2 };

static inline Lut3D get_neg_lut(int choice){
    using namespace EmbeddedLUTs;
    switch(choice){
        default:
        case Neg_Cthulhu:    return {LUT_NEG_CTHULHU, LUT_SIZE};
        case Neg_Lilith:     return {LUT_NEG_LILITH, LUT_SIZE};
        case Neg_Tsathoggua: return {LUT_NEG_TSATHOGGUA, LUT_SIZE};
        case Neg_Yig:        return {LUT_NEG_YIG, LUT_SIZE};
    }
}
static inline Lut3D get_sep_lut(int choice){
    using namespace EmbeddedLUTs;
    switch(choice){
        default:
        case Sep_Hydra: return {LUT_SEP_HYDRA, LUT_SIZE};
        case Sep_Oorn:  return {LUT_SEP_OORN, LUT_SIZE};
        case Sep_Zhar:  return {LUT_SEP_ZHAR, LUT_SIZE};
    }
}
static inline Lut3D get_print_lut(int /*choice*/){
    using namespace EmbeddedLUTs;
    return {LUT_PRINT_KODAK, LUT_SIZE};
}

class FilmPipelineProcessor : public OFX::ImageProcessor {
public:
    FilmPipelineProcessor(OFX::ImageEffect& instance)
    : OFX::ImageProcessor(instance)
    , _srcImg(nullptr), _dstImg(nullptr)
    {}

    void setSrcImg(const OFX::Image* img){ _srcImg = img; }
    void setDstImg(OFX::Image* img){ _dstImg = img; }

    // Params (filled per render)
    int inGamut = 15; // DaVinciWG
    int inOetf  = 1;  // DaVinci Intermediate

    bool negEnable = true;
    int  negChoice = 0;
    float negBlend = 0.8f;

    bool sepEnable = true;
    int  sepChoice = 0;
    float sepBlend = 0.5f;

    bool printEnable = true;
    int  printChoice = 0; // Kodak only
    float printBlend = 0.5f;

private:
    const OFX::Image* _srcImg;
    OFX::Image* _dstImg;

    void multiThreadProcessImages(OfxRectI procWindow) override {
        if(!_srcImg || !_dstImg) return;

        const int nComp = _dstImg->getPixelComponentCount();
        const OFX::BitDepthEnum bd = _dstImg->getPixelDepth();
        if(bd != OFX::eBitDepthFloat || (nComp != 4 && nComp != 3)) {
            OFX::throwSuiteStatusException(kOfxStatErrUnsupported);
        }

        const Lut3D negLut   = get_neg_lut(negChoice);
        const Lut3D sepLut   = get_sep_lut(sepChoice);
        const Lut3D printLut = get_print_lut(printChoice);

        for(int y = procWindow.y1; y < procWindow.y2; ++y) {
            if(_effect.abort()) break;

            float* dstPix = (float*)_dstImg->getPixelAddress(procWindow.x1, y);
            for(int x = procWindow.x1; x < procWindow.x2; ++x) {
                const float* srcPix = (const float*)_srcImg->getPixelAddress(x, y);

                float3 rgb_in = { srcPix[0], srcPix[1], srcPix[2] };
                float  a_in   = (nComp == 4) ? srcPix[3] : 1.0f;

                // 1) Color-manage to DWG+DI
                float3 dwg_di = input_to_dwg_intermediate(rgb_in, inGamut, inOetf);

                // 2) Negative (luma LUT) in DWG+DI
                float3 after_neg = dwg_di;
                if(negEnable){
                    float Y = luma_rec709(dwg_di);
                    float3 neutral = {Y, Y, Y};
                    float3 lutOut = lut_sample_tetra(negLut, neutral);
                    float Y2 = luma_rec709(lutOut);
                    float scale = (Y > 1e-6f) ? (Y2 / Y) : 1.0f;
                    float3 scaled = dwg_di * scale;
                    after_neg = lerp3(dwg_di, scaled, clampf(negBlend, 0.0f, 1.0f));
                }

                // 3) Color separation LUT in DWG+DI
                float3 after_sep = after_neg;
                if(sepEnable){
                    float3 sepOut = lut_sample_tetra(sepLut, after_neg);
                    after_sep = lerp3(after_neg, sepOut, clampf(sepBlend, 0.0f, 1.0f));
                }

                // 4) Print stage in output space (Rec709 2.4 baseline vs Kodak LUT)
                float3 baseline = dwg_di_to_rec709_24(after_sep);

                float3 out_rgb = baseline;
                if(printEnable){
                    // Kodak LUT assumed to take DWG+DI and output Rec709-ish
                    float3 kodak = lut_sample_tetra(printLut, after_sep);
                    out_rgb = lerp3(baseline, kodak, clampf(printBlend, 0.0f, 1.0f));
                }

                dstPix[0] = out_rgb.x;
                dstPix[1] = out_rgb.y;
                dstPix[2] = out_rgb.z;
                if(nComp == 4) dstPix[3] = a_in;

                dstPix += nComp;
            }
        }
    }
};

class FilmPipelineEffect : public OFX::ImageEffect {
public:
    FilmPipelineEffect(OfxImageEffectHandle handle)
    : ImageEffect(handle)
    , _srcClip(nullptr), _dstClip(nullptr)
    , _pInGamut(nullptr), _pInOetf(nullptr)
    , _pNegEnable(nullptr), _pNegLut(nullptr), _pNegBlend(nullptr)
    , _pSepEnable(nullptr), _pSepStyle(nullptr), _pSepBlend(nullptr)
    , _pPrintEnable(nullptr), _pPrintLut(nullptr), _pPrintBlend(nullptr)
    {
        _dstClip = fetchClip(kOfxImageEffectOutputClipName);
        _srcClip = fetchClip(kOfxImageEffectSimpleSourceClipName);

        _pInGamut = fetchChoiceParam(kParamInGamut);
        _pInOetf  = fetchChoiceParam(kParamInOETF);

        _pNegEnable = fetchBooleanParam(kParamNegEnable);
        _pNegLut    = fetchChoiceParam(kParamNegLut);
        _pNegBlend  = fetchDoubleParam(kParamNegBlend);

        _pSepEnable = fetchBooleanParam(kParamSepEnable);
        _pSepStyle  = fetchChoiceParam(kParamSepStyle);
        _pSepBlend  = fetchDoubleParam(kParamSepBlend);

        _pPrintEnable = fetchBooleanParam(kParamPrintEnable);
        _pPrintLut    = fetchChoiceParam(kParamPrintLut);
        _pPrintBlend  = fetchDoubleParam(kParamPrintBlend);
    }

private:
    OFX::Clip* _srcClip;
    OFX::Clip* _dstClip;

    OFX::ChoiceParam* _pInGamut;
    OFX::ChoiceParam* _pInOetf;

    OFX::BooleanParam* _pNegEnable;
    OFX::ChoiceParam*  _pNegLut;
    OFX::DoubleParam*  _pNegBlend;

    OFX::BooleanParam* _pSepEnable;
    OFX::ChoiceParam*  _pSepStyle;
    OFX::DoubleParam*  _pSepBlend;

    OFX::BooleanParam* _pPrintEnable;
    OFX::ChoiceParam*  _pPrintLut;
    OFX::DoubleParam*  _pPrintBlend;

    void render(const OFX::RenderArguments& args) override {
        std::unique_ptr<OFX::Image> dst(_dstClip->fetchImage(args.time));
        std::unique_ptr<const OFX::Image> src(_srcClip->fetchImage(args.time));

        if(!dst || !src) OFX::throwSuiteStatusException(kOfxStatFailed);

        FilmPipelineProcessor proc(*this);
        proc.setDstImg(dst.get());
        proc.setSrcImg(src.get());
        proc.setRenderWindow(args.renderWindow);

        // Read params
                _pInGamut->getValue(proc.inGamut);
        _pInOetf->getValue(proc.inOetf);
proc.negEnable = _pNegEnable->getValue();
                _pNegLut->getValue(proc.negChoice);
proc.negBlend  = (float)_pNegBlend->getValue();

        proc.sepEnable = _pSepEnable->getValue();
                _pSepStyle->getValue(proc.sepChoice);
proc.sepBlend  = (float)_pSepBlend->getValue();

        proc.printEnable = _pPrintEnable->getValue();
                _pPrintLut->getValue(proc.printChoice);
proc.printBlend  = (float)_pPrintBlend->getValue();

        proc.process();
    }
};

class FilmPipelineFactory : public OFX::PluginFactoryHelper<FilmPipelineFactory> {
public:
    FilmPipelineFactory(): OFX::PluginFactoryHelper<FilmPipelineFactory>(kPluginIdentifier, kPluginVersionMajor, kPluginVersionMinor) {}
    void describe(OFX::ImageEffectDescriptor& desc) override {
        desc.setLabels(kPluginName, kPluginName, kPluginName);
        desc.setPluginGrouping(kPluginGrouping);
        desc.setPluginDescription(kPluginDescription);

        desc.addSupportedContext(OFX::eContextFilter);
        desc.addSupportedBitDepth(OFX::eBitDepthFloat);

        desc.setSingleInstance(false);
        desc.setHostFrameThreading(false);
        desc.setSupportsMultiResolution(true);
        desc.setSupportsTiles(true);
        desc.setRenderTwiceAlways(false);
        desc.setSupportsMultipleClipPARs(false);

        // Clips
        OFX::ClipDescriptor* srcClip = desc.defineClip(kOfxImageEffectSimpleSourceClipName);
        srcClip->addSupportedComponent(OFX::ePixelComponentRGBA);
        srcClip->addSupportedComponent(OFX::ePixelComponentRGB);
        srcClip->setTemporalClipAccess(false);
        srcClip->setSupportsTiles(true);

        OFX::ClipDescriptor* dstClip = desc.defineClip(kOfxImageEffectOutputClipName);
        dstClip->addSupportedComponent(OFX::ePixelComponentRGBA);
        dstClip->addSupportedComponent(OFX::ePixelComponentRGB);
        dstClip->setSupportsTiles(true);
    }

    void describeInContext(OFX::ImageEffectDescriptor& desc, OFX::ContextEnum) override {
        OFX::PageParamDescriptor* page = desc.definePageParam("Controls");

        // Input
        {
            OFX::ChoiceParamDescriptor* p = desc.defineChoiceParam(kParamInGamut);
            p->setLabel("Input Gamut");
            p->appendOption("XYZ");
            p->appendOption("ACES 2065-1 (AP0)");
            p->appendOption("ACEScg (AP1)");
            p->appendOption("P3 D65");
            p->appendOption("Rec.2020");
            p->appendOption("Rec.709");
            p->appendOption("Arri Wide Gamut 3");
            p->appendOption("Arri Wide Gamut 4");
            p->appendOption("RED Wide Gamut RGB");
            p->appendOption("Sony SGamut3");
            p->appendOption("Sony SGamut3Cine");
            p->appendOption("Panasonic V-Gamut");
            p->appendOption("Blackmagic Wide Gamut");
            p->appendOption("Filmlight E-Gamut");
            p->appendOption("Filmlight E-Gamut2");
            p->appendOption("DaVinci Wide Gamut");
            p->setDefault(15);
            if(page) page->addChild(*p);

            OFX::ChoiceParamDescriptor* t = desc.defineChoiceParam(kParamInOETF);
            t->setLabel("Input Transfer");
            t->appendOption("Linear");
            t->appendOption("DaVinci Intermediate");
            t->appendOption("Filmlight T-Log");
            t->appendOption("ACEScct");
            t->appendOption("Arri LogC3");
            t->appendOption("Arri LogC4");
            t->appendOption("RedLog3G10");
            t->appendOption("Panasonic V-Log");
            t->appendOption("Sony S-Log3");
            t->appendOption("Fuji F-Log2");
            t->setDefault(1);
            if(page) page->addChild(*t);
        }

        // Negative
        {
            OFX::BooleanParamDescriptor* e = desc.defineBooleanParam(kParamNegEnable);
            e->setLabel("Negative Enable");
            e->setDefault(true);
            if(page) page->addChild(*e);

            OFX::ChoiceParamDescriptor* p = desc.defineChoiceParam(kParamNegLut);
            p->setLabel("Negative LUT");
            p->appendOption("Cthulhu");
            p->appendOption("Lilith");
            p->appendOption("Tsathoggua");
            p->appendOption("Yig");
            p->setDefault(0);
            if(page) page->addChild(*p);

            OFX::DoubleParamDescriptor* b = desc.defineDoubleParam(kParamNegBlend);
            b->setLabel("Negative Blend");
            b->setDefault(0.8);
            b->setRange(0.0, 1.0);
            b->setDisplayRange(0.0, 1.0);
            if(page) page->addChild(*b);
        }

        // Color Separation
        {
            OFX::BooleanParamDescriptor* e = desc.defineBooleanParam(kParamSepEnable);
            e->setLabel("Color Separation Enable");
            e->setDefault(true);
            if(page) page->addChild(*e);

            OFX::ChoiceParamDescriptor* p = desc.defineChoiceParam(kParamSepStyle);
            p->setLabel("Color Separation Style");
            p->appendOption("Hydra");
            p->appendOption("Oorn");
            p->appendOption("Zhar");
            p->setDefault(0);
            if(page) page->addChild(*p);

            OFX::DoubleParamDescriptor* b = desc.defineDoubleParam(kParamSepBlend);
            b->setLabel("Color Separation Blend");
            b->setDefault(0.5);
            b->setRange(0.0, 1.0);
            b->setDisplayRange(0.0, 1.0);
            if(page) page->addChild(*b);
        }

        // Print
        {
            OFX::BooleanParamDescriptor* e = desc.defineBooleanParam(kParamPrintEnable);
            e->setLabel("Print Enable");
            e->setDefault(true);
            if(page) page->addChild(*e);

            OFX::ChoiceParamDescriptor* p = desc.defineChoiceParam(kParamPrintLut);
            p->setLabel("Print LUT");
            p->appendOption("Kodak");
            p->setDefault(0);
            if(page) page->addChild(*p);

            OFX::DoubleParamDescriptor* b = desc.defineDoubleParam(kParamPrintBlend);
            b->setLabel("Print Blend");
            b->setDefault(0.5);
            b->setRange(0.0, 1.0);
            b->setDisplayRange(0.0, 1.0);
            if(page) page->addChild(*b);
        }
    }

    OFX::ImageEffect* createInstance(OfxImageEffectHandle handle, OFX::ContextEnum) override {
        return new FilmPipelineEffect(handle);
    }
};

static FilmPipelineFactory pFactory;

void OFX::Plugin::getPluginIDs(OFX::PluginFactoryArray& ids) {
    ids.push_back(&pFactory);
}
