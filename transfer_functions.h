#pragma once
#include <cmath>
#include <algorithm>

static inline float _log2f(float x){ return std::log(x)/std::log(2.0f); }
static inline float _exp2f(float x){ return std::exp(x*std::log(2.0f)); }
static inline float _expf(float x){ return std::exp(x); }
static inline float _powf(float a, float b){ return std::pow(a,b); }
static inline float _exp10f(float x){ return std::pow(10.0f, x); }
static inline float _fmaxf(float a,float b){ return std::max(a,b); }
static inline float _fminf(float a,float b){ return std::min(a,b); }


// Linearization functions (encoded -> linear), adapted from OpenDRT
/* OETF Linearization Transfer Functions ---------------------------------------- */

static inline float oetf_davinci_intermediate(float x) {
    return x <= 0.02740668f ? x/10.44426855f : _exp2f(x/0.07329248f - 7.0f) - 0.0075f;
}
static inline float oetf_filmlight_tlog(float x) {
  return x < 0.075f ? (x-0.075f)/16.184376489665897f : _expf((x - 0.5520126568606655f)/0.09232902596577353f) - 0.0057048244042473785f;
}
static inline float oetf_acescct(float x) {
  return x <= 0.155251141552511f ? (x - 0.0729055341958355f)/10.5402377416545f : _exp2f(x*17.52f - 9.72f);
}
static inline float oetf_arri_logc3(float x) {
  return x < 5.367655f*0.010591f + 0.092809f ? (x - 0.092809f)/5.367655f : (_exp10f((x - 0.385537f)/0.247190f) - 0.052272f)/5.555556f;
}
static inline float oetf_arri_logc4(float x) {
  return x < -0.7774983977293537f ? x*0.3033266726886969f - 0.7774983977293537f : (_exp2f(14.0f*(x - 0.09286412512218964f)/0.9071358748778103f + 6.0f) - 64.0f)/2231.8263090676883f;
}
static inline float oetf_red_log3g10(float x) {
  return x < 0.0f ? (x/15.1927f) - 0.01f : (_exp10f(x/0.224282f) - 1.0f)/155.975327f - 0.01f;
}
static inline float oetf_panasonic_vlog(float x) {
  return x < 0.181f ? (x - 0.125f)/5.6f : _exp10f((x - 0.598206f)/0.241514f) - 0.00873f;
}
static inline float oetf_sony_slog3(float x) {
  return x < 171.2102946929f/1023.0f ? (x*1023.0f - 95.0f)*0.01125f/(171.2102946929f - 95.0f) : (_exp10f(((x*1023.0f - 420.0f)/261.5f))*(0.18f + 0.01f) - 0.01f);
}
static inline float oetf_fujifilm_flog2(float x) {
  return x < 0.100686685370811f ? (x - 0.092864f)/8.799461f : (_exp10f(((x - 0.384316f)/0.245281f))/5.555556f - 0.064829f/5.555556f);
}


static inline float decode_input_oetf(int idx, float x){
  switch(idx){
    case 0: return x;
    case 1: return oetf_davinci_intermediate(x);
    case 2: return oetf_filmlight_tlog(x);
    case 3: return oetf_acescct(x);
    case 4: return oetf_arri_logc3(x);
    case 5: return oetf_arri_logc4(x);
    case 6: return oetf_red_log3g10(x);
    case 7: return oetf_panasonic_vlog(x);
    case 8: return oetf_sony_slog3(x);
    case 9: return oetf_fujifilm_flog2(x);
    default: return x;
  }
}

static inline float encode_davinci_intermediate(float y){
  if(y <= 0.002624088021941948f) return y * 10.44426855f;
  return 0.07329248f * (_log2f(std::max(y + 0.0075f, 1e-12f)) + 7.0f);
}

static inline float decode_davinci_intermediate(float x){ return oetf_davinci_intermediate(x); }
static inline float encode_rec709_24(float x){ x = std::max(x, 0.0f); return std::pow(x, 1.0f/2.4f); }
