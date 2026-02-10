#pragma once
#include <array>
#include <algorithm>

struct Mat3 {
  float m[3][3];
  Mat3() {
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        m[r][c] = 0.0f;
      }
    }
  }
  Mat3(float a00, float a01, float a02,
       float a10, float a11, float a12,
       float a20, float a21, float a22) {
    m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
    m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
    m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
  }
};


static inline Mat3 mat_identity(){ return Mat3(1.0f,0.0f,0.0f, 0.0f,1.0f,0.0f, 0.0f,0.0f,1.0f); }

static inline float det3(const Mat3& A){
  return A.m[0][0]*(A.m[1][1]*A.m[2][2]-A.m[1][2]*A.m[2][1])
       - A.m[0][1]*(A.m[1][0]*A.m[2][2]-A.m[1][2]*A.m[2][0])
       + A.m[0][2]*(A.m[1][0]*A.m[2][1]-A.m[1][1]*A.m[2][0]);
}
static inline Mat3 inv3(const Mat3& A){
  float d = det3(A);
  if(d==0.0f) return mat_identity();
  float id = 1.0f/d;
  Mat3 R;
  R.m[0][0] =  (A.m[1][1]*A.m[2][2]-A.m[1][2]*A.m[2][1])*id;
  R.m[0][1] = -(A.m[0][1]*A.m[2][2]-A.m[0][2]*A.m[2][1])*id;
  R.m[0][2] =  (A.m[0][1]*A.m[1][2]-A.m[0][2]*A.m[1][1])*id;
  R.m[1][0] = -(A.m[1][0]*A.m[2][2]-A.m[1][2]*A.m[2][0])*id;
  R.m[1][1] =  (A.m[0][0]*A.m[2][2]-A.m[0][2]*A.m[2][0])*id;
  R.m[1][2] = -(A.m[0][0]*A.m[1][2]-A.m[0][2]*A.m[1][0])*id;
  R.m[2][0] =  (A.m[1][0]*A.m[2][1]-A.m[1][1]*A.m[2][0])*id;
  R.m[2][1] = -(A.m[0][0]*A.m[2][1]-A.m[0][1]*A.m[2][0])*id;
  R.m[2][2] =  (A.m[0][0]*A.m[1][1]-A.m[0][1]*A.m[1][0])*id;
  return R;
}
static inline std::array<float,3> mat_vec(const Mat3& A, const std::array<float,3>& v){
  return {A.m[0][0]*v[0]+A.m[0][1]*v[1]+A.m[0][2]*v[2],
          A.m[1][0]*v[0]+A.m[1][1]*v[1]+A.m[1][2]*v[2],
          A.m[2][0]*v[0]+A.m[2][1]*v[1]+A.m[2][2]*v[2]};
}

static const Mat3 matrix_ap0_to_xyz = Mat3(0.938630949f, -0.00574192055f, 0.0175668989f, 0.338093595f, 0.727213903f, -0.0653074977f, 0.000723121511f, 0.000818441849f, 1.08751619f);
static const Mat3 matrix_ap1_to_xyz = Mat3(0.652418718f, 0.127179926f, 0.170857284f, 0.268064059f, 0.672464479f, 0.0594714618f, -0.00546992851f, 0.00518279998f, 1.08934488f);
static const Mat3 matrix_p3d65_to_xyz = Mat3(0.486570949f, 0.265667693f, 0.198217285f, 0.228974564f, 0.691738522f, 0.0792869141f, 0.0f, 0.0451133819f, 1.04394437f);
static const Mat3 matrix_rec2020_to_xyz = Mat3(0.636958048f, 0.144616904f, 0.168880975f, 0.262700212f, 0.677998072f, 0.0593017165f, 0.0f, 0.028072693f, 1.06098506f);
static const Mat3 matrix_rec709_to_xyz = Mat3(0.412390799f, 0.357584339f, 0.180480788f, 0.212639006f, 0.715168679f, 0.0721923154f, 0.0193308187f, 0.11919478f, 0.950532152f);
static const Mat3 matrix_arriwg3_to_xyz = Mat3(0.638007619f, 0.214703856f, 0.0977444514f, 0.291953779f, 0.823841042f, -0.115794821f, 0.00279827903f, -0.0670342357f, 1.15329371f);
static const Mat3 matrix_arriwg4_to_xyz = Mat3(0.70485832f, 0.129760295f, 0.115837311f, 0.254524176f, 0.781477733f, -0.0360019091f, 0.0f, 0.0f, 1.08905775f);
static const Mat3 matrix_redwg_to_xyz = Mat3(0.735275246f, 0.0686094106f, 0.146571271f, 0.286694099f, 0.842979134f, -0.129673234f, -0.0796808569f, -0.347343217f, 1.51608182f);
static const Mat3 matrix_sonysgamut3_to_xyz = Mat3(0.706482713f, 0.12880105f, 0.115172164f, 0.270979671f, 0.786606411f, -0.057586082f, -0.00967784539f, 0.00460003749f, 1.09413556f);
static const Mat3 matrix_sonysgamut3cine_to_xyz = Mat3(0.599083921f, 0.248925516f, 0.10244649f, 0.21507582f, 0.885068502f, -0.100144322f, -0.0320658495f, -0.0276583907f, 1.14878199f);
static const Mat3 matrix_vgamut_to_xyz = Mat3(0.67964447f, 0.152211412f, 0.118600045f, 0.26068555f, 0.774894463f, -0.0355800134f, -0.00931019822f, -0.00461246704f, 1.10298042f);
static const Mat3 matrix_bmdwg_to_xyz = Mat3(0.606538368f, 0.220412735f, 0.123504823f, 0.26799294f, 0.832748409f, -0.100741349f, -0.0294425542f, -0.0866124303f, 1.20511274f);
static const Mat3 matrix_egamut_to_xyz = Mat3(0.70539685f, 0.164041328f, 0.0810177487f, 0.280130724f, 0.820206642f, -0.100337366f, -0.103781512f, -0.072907257f, 1.26574652f);
static const Mat3 matrix_egamut2_to_xyz = Mat3(0.7364777f, 0.130739651f, 0.0832385758f, 0.275069984f, 0.82801779f, -0.103087775f, -0.124225154f, -0.0871597674f, 1.30044267f);
static const Mat3 matrix_davinciwg_to_xyz = Mat3(0.700622392f, 0.148774815f, 0.10105872f, 0.274118511f, 0.873631896f, -0.147750407f, -0.0989629129f, -0.137895325f, 1.32591599f);

static inline Mat3 inputGamutToXYZ(int idx){
  switch(idx){
    case 0: return mat_identity(); // XYZ
    case 1: return matrix_ap0_to_xyz; // ACES 2065-1
    case 2: return matrix_ap1_to_xyz; // ACEScg
    case 3: return matrix_p3d65_to_xyz; // P3D65
    case 4: return matrix_rec2020_to_xyz; // Rec2020
    case 5: return matrix_rec709_to_xyz; // Rec709
    case 6: return matrix_arriwg3_to_xyz; // ArriWG3
    case 7: return matrix_arriwg4_to_xyz; // ArriWG4
    case 8: return matrix_redwg_to_xyz; // RedWG
    case 9: return matrix_sonysgamut3_to_xyz; // SonySGamut3
    case 10: return matrix_sonysgamut3cine_to_xyz; // SonySGamut3Cine
    case 11: return matrix_vgamut_to_xyz; // PanasonicVGamut
    case 12: return matrix_bmdwg_to_xyz; // BMDWG
    case 13: return matrix_egamut_to_xyz; // FilmlightEGamut
    case 14: return matrix_egamut2_to_xyz; // FilmlightEGamut2
    case 15: return matrix_davinciwg_to_xyz; // DaVinciWG
    default: return mat_identity();
  }
}

static inline Mat3 xyzToDaVinciWG(){ static Mat3 inv = inv3(matrix_davinciwg_to_xyz); return inv; }
static inline Mat3 xyzToRec709(){ static Mat3 inv = inv3(matrix_rec709_to_xyz); return inv; }


