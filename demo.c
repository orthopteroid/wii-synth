// 2014 orthopteroid@gmail.com
// kudos to Youth Uprising, Perry Cook and the SDL gang
// compile with:
// gcc demo.c -o demo -ffast-math -fomit-frame-pointer -msseregparm -mfpmath=sse -msse3 `sdl-config --cflags --libs` -lwiiuse -lm

#include <stdio.h>
#include <math.h>

#include <SDL/SDL.h>
#include "wiiuse.h"                     /* for wiimote_t, classic_ctrl_t, etc */

#include <emmintrin.h>

#include <time.h>

#define MAX_WIIMOTES 1

//#define ENABLE_RENDERTIMING
//#define ENABLE_RENDERCHECKS
//#define ENABLE_RENDERDUMP
//#define ENABLE_RENDERDUMP_STOP

#define SPINLOCK_SET(sl)	do{ (sl) = 1; } while(0)
#define SPINLOCK_CLEAR(sl)	do{ (sl) = 0; } while(0)
#define SPINLOCK_WAIT(sl)	do{ while( sl ) { ; } } while(0)

//////////////////////////

static int app_stop = 0;

//////////////////////////

static wiimote** wii_wiimotes = 0;
static int wii_connected = 0;

//////////////////////////

#define BUF_TYPE_SDL	AUDIO_S16SYS
#define BUF_TYPE		Uint16
#define BUF_SCALE		(1 << (sizeof(BUF_TYPE) * 8 - 1) - 1)

void buf_mixerCallback(void *userdata, Uint8 *stream, int len); // fwd decl

SDL_AudioSpec sdl_configrequest =
{
	22050 /*freq*/, BUF_TYPE_SDL /*format*/, 1 /*channels*/, 0 /*silence*/,
	4096 /*samples*/, 0 /*padding*/, 4096 /*size*/, buf_mixerCallback, 0 /*userdata*/
};
SDL_AudioSpec sdl_config;

#define BUF_COUNT 2 // power of 2, increase if skipping occurs
#define BUF_WAIT_FIRSTUSE()			do { while( !buf_writeID ) { ; } } while(0)
#define BUF_WAIT_READCOMPLETE()		do { while( buf_playID == buf_writeID ) { ; } } while(0)
#define BUF_INC_WRITE()				do { buf_writeID++; } while(0)
#define BUF_INC_READ()				do { buf_playID++; } while(0)
#define BUF_GET_WRITE()				( buf_writeID & (BUF_COUNT - 1) )
#define BUF_GET_PLAY()				( buf_playID  & (BUF_COUNT - 1) )

static BUF_TYPE *sdl_buffer[ BUF_COUNT ];
static Uint32 buf_playID = BUF_COUNT - 1;
static Uint32 buf_writeID = 0;

void buf_init()
{
	int i = 0;
	for(; i < BUF_COUNT; i++ )
	{
		sdl_buffer[i] = (BUF_TYPE*)malloc( sdl_config.size );
		memset( (void*)sdl_buffer[i], 0, sdl_config.size );
	}
}

static Uint32 buf_transferCount = 0;
static Uint32 buf_skipCount = 0;
static Uint32 buf_chopCount = 0;

void buf_cleanup()
{
	int i = 0;
	for(; i < BUF_COUNT; i++ )
	{
		free( (void*)sdl_buffer[i] );
	}

	printf("%u skips (render took too long)\n", buf_skipCount );
	printf("%u chops (partial-buffer transfers)\n", buf_chopCount );
}

void buf_mixerCallback(void *userdata, Uint8 *stream, int len)
{
	buf_transferCount++;
	if( len != sdl_config.size ) { buf_chopCount++; }
	if( buf_playID - BUF_COUNT > buf_writeID ) { buf_skipCount++; }

	SDL_MixAudio( stream, (Uint8*)sdl_buffer[ BUF_GET_PLAY() ], len, SDL_MIX_MAXVOLUME );
	BUF_INC_READ();
}

//////////////////////////

// http://stackoverflow.com/questions/3388134/rdtsc-accuracy-across-cpu-cores
unsigned long perf_tsc()
{
	unsigned long v;
	__asm__ volatile (".byte 0x0f, 0x31; orl $1, %%eax; movl %%eax, %0;" : "=r" (v) : : "%eax" /*clobbered*/ ); // lsb of rdtsc (in eax) made odd and stored
	return v;
}

//////////////////////////

static float synth_filterNetwork[100];

static SDL_Thread *synth_threadID = 0;
static int synth_threadCode = 0;
static int synth_renderSpinlock = 0;

void synth_renderEvent( Sint8 c, float rollIn, float pitchIn, float yawIn ); // fwd decl
float synth_render( Uint32 uT ); // fwd decl

static unsigned long synth_loopCycles = 0;
static unsigned long synth_loopIter = 0;

int synth_renderThread(void *ptr)
{
	Uint32 t = 0, i = 0;
	Uint16 err_nan, err_range;

	while( !app_stop )
	{
		err_nan = err_range = 0;

		unsigned long loop_start = perf_tsc();

		SPINLOCK_SET( synth_renderSpinlock );
		for( i=0; i<sdl_config.samples; i++, t++ )
		{
			float Q = synth_render( t );

#if defined(ENABLE_RENDERCHECKS)

			if( isnan( Q ) ) err_nan = 1;
			else if( !((-1.f <= Q) & (Q <= 1.f)) ) err_range = 1;

#endif // ENABLE_RENDERCHECKS

			BUF_TYPE iQ = (BUF_TYPE)( (float)BUF_SCALE + (float)BUF_SCALE * Q );
			sdl_buffer[ BUF_GET_WRITE() ][ i ] = iQ;

#if defined(ENABLE_RENDERDUMP)

			printf("%4u ",iQ);
			static int c = 0;
			if( ++c > 21 ) { c=0; putchar('\n');}

#if defined(ENABLE_RENDERDUMP_STOP)

			static int z = 0; if( ++z > 50000 ) { app_stop = 1; }

#endif // ENABLE_RENDERDUMP_STOP

#endif // ENABLE_RENDERDUMP

		}
		SPINLOCK_CLEAR( synth_renderSpinlock );

		synth_loopCycles += perf_tsc() - loop_start; synth_loopIter++;

#if defined(ENABLE_RENDERCHECKS)

		if( err_nan )			printf("N\n"); // nan error
		else if( err_range )	printf("R\n"); // range error

#endif // ENABLE_RENDERCHECKS

		BUF_INC_WRITE();
		BUF_WAIT_READCOMPLETE();
	}

    return 0;
}

//////////////////////////

#define FLT_NEW(flt) \
	const Uint32 (flt) = (__COUNTER__)

#define FLT_OUTPUT(flt) (synth_filterNetwork[flt])

#define FLT_SET(flt, value) \
	do { \
		synth_filterNetwork[flt] = (value); \
	} while(0)

#define FLT_AMP(flt, incoef, in, biascoef) \
	do { \
		synth_filterNetwork[flt] = (incoef) * synth_filterNetwork[in] + (biascoef); \
	} while(0)

#define FLT_MIX2(flt, in0, in1) \
	do { \
		synth_filterNetwork[flt]  = synth_filterNetwork[in0]; \
		synth_filterNetwork[flt] += synth_filterNetwork[in1]; \
		synth_filterNetwork[flt] /= 2.f; \
	} while(0)

#define FLT_MIX3(flt, in0, in1, in2) \
	do { \
		synth_filterNetwork[flt]  = synth_filterNetwork[in0]; \
		synth_filterNetwork[flt] += synth_filterNetwork[in1]; \
		synth_filterNetwork[flt] += synth_filterNetwork[in2]; \
		synth_filterNetwork[flt] /= 3.f; \
	} while(0)

#define FLT_MIX2_CROSS(flt, coefN1toP1, in0, in1) \
	do { \
		float coef0to1 = coefN1toP1 * .5f + .5f; /* convert [-1,+1] to [0,1] */ \
		synth_filterNetwork[flt] = (coef0to1) * synth_filterNetwork[in0] + (1.f - (coef0to1)) * synth_filterNetwork[in1]; \
	} while(0)

/////////////////////////////

#define FLT_SWEEP_NEW(flt) \
	const Uint32 flt = (__COUNTER__), flt ## _storage = (__COUNTER__ + __COUNTER__ + __COUNTER__ + __COUNTER__)

#define FLT_SWEEP_WRAP_CONFIG(flt, step, minval, maxval) \
	do { \
		synth_filterNetwork[flt+1] = 0.f; \
		synth_filterNetwork[flt+2] = step; \
		synth_filterNetwork[flt+3] = (minval); \
		synth_filterNetwork[flt+4] = (maxval); \
	} while(0)

#define FLT_SWEEP_WRAP(flt, rate) \
	math_interp_wrap( &(synth_filterNetwork[flt]), (rate), synth_filterNetwork[flt+3], synth_filterNetwork[flt+4] )

#define FLT_SWEEP_PONG_CONFIG(flt, step, minval, maxval) \
	do { \
		synth_filterNetwork[flt+1] = math_sign(step); \
		synth_filterNetwork[flt+2] = math_fabs(step); \
		synth_filterNetwork[flt+3] = (minval); \
		synth_filterNetwork[flt+4] = (maxval); \
	} while(0)

#define FLT_SWEEP_PONG(flt) \
	math_interp_pong( &(synth_filterNetwork[flt]), &(synth_filterNetwork[flt+1]), synth_filterNetwork[flt+2], synth_filterNetwork[flt+3], synth_filterNetwork[flt+4] )

#define FLT_SWEEP_PONG_RATED(flt, rate) \
	math_interp_pong( &(synth_filterNetwork[flt]), &(synth_filterNetwork[flt+1]), (rate), synth_filterNetwork[flt+3], synth_filterNetwork[flt+4] )

/////////////////////////////

#define FLT_INTERP_NEW(flt) \
	const Uint32 flt = (__COUNTER__), flt ## _storage = (__COUNTER__ + __COUNTER__)

#define FLT_INTERP_TRIGGER(flt, velocity, target) \
	do { \
		synth_filterNetwork[flt+1] = (velocity) * ( (target) > synth_filterNetwork[flt] ? +1.f : -1.f ); \
		synth_filterNetwork[flt+2] = (target); \
	} while(0)

#define FLT_INTERP(flt) \
	math_interp_clamp( &(synth_filterNetwork[flt]), synth_filterNetwork[flt+1], synth_filterNetwork[flt+2] )

/////////////////////////////

#define FLT_ALLPASS(flt, in, coef0to1) \
	do { \
		float apf_in0 = synth_filterNetwork[in][0]; \
		float apf_in1 = synth_filterNetwork[flt][1]; \
		float apf_out1 = (1.f - (coef0to1)) * apf_in1 + (coef0to1) * apf_in0; \
		float apf_out0 = (1.f - (coef0to1)) * apf_in0 + (coef0to1) * apf_out1; \
		synth_filterNetwork[flt][0] = apf_out0; \
		synth_filterNetwork[flt][1] = apf_out1; /* storage */ \
	} while(0)

/////////////////////////////

#define FLT_IIRF_NEW(flt) \
	const Uint32 flt = (__COUNTER__), flt ## _storage = (__COUNTER__ + __COUNTER__ + __COUNTER__ + __COUNTER__)

#define FLT_IIRRF_TRIGGER(flt, samplerate, kcoef, mcoef, rcoef) \
	do { \
		float samples = 50.f; \
		float oT = 1.f / ((samplerate) * samples); \
		float denom = (mcoef) + oT * (rcoef) + oT * oT * (kcoef); \
		synth_filterNetwork[flt] = 0.f; \
		synth_filterNetwork[flt+1] = 1.f; \
		synth_filterNetwork[flt+2] = 0.f; \
		synth_filterNetwork[flt+3] = (2.f * (mcoef) + oT * (rcoef)) / denom; \
		synth_filterNetwork[flt+4] = -(mcoef) / denom; \
	} while(0)

#define FLT_IIRRF(flt) \
	do { \
		Uint16 i; \
		for( i=0; i<50; i++ ) \
		{ \
			float tmp = synth_filterNetwork[flt+3] * synth_filterNetwork[flt+1] + synth_filterNetwork[flt+4] * synth_filterNetwork[flt+2]; \
			synth_filterNetwork[flt+2] = synth_filterNetwork[flt+1]; \
			synth_filterNetwork[flt+1] = tmp; \
			synth_filterNetwork[flt] = synth_filterNetwork[flt+1] - synth_filterNetwork[flt+2]; \
		} \
	} while(0)

//////////////////////////////

#define FLT_FUNC_NEW(flt) \
	const Uint32 flt = (__COUNTER__)

#define FLT_FUNC_WHITE(flt, amp) \
	do { \
		synth_filterNetwork[flt] = (amp) * ( 2.f * math_rand() - 1.f ); \
	} while(0)

#define FLT_FUNC_SIN(flt, amp, vlf_osc3_angle) \
	do { \
		synth_filterNetwork[flt] = (amp) * math_sin( vlf_osc3_angle ); \
	} while(0)

#define FLT_FUNC_SQUARE(flt, amp, vlf_osc3_angle) \
	do { \
		synth_filterNetwork[flt] = (amp) * math_sign( vlf_osc3_angle ); \
	} while(0)

#define FLT_FUNC_SAW(flt, amp, vlf_osc3_angle) \
	do { \
		synth_filterNetwork[flt] = (amp) * math_saw( vlf_osc3_angle ); \
	} while(0)

#define FLT_FUNC_SIMILAR(flt, in0, in1) \
	do { \
		synth_filterNetwork[flt] = math_percent_similar( synth_filterNetwork[in0], synth_filterNetwork[in1] ); \
	} while(0)

//////////////////////////

// the optimizations don't seem to work with a multithreaded sdl app, so we turn them on locally only
#pragma GCC push_options
#pragma GCC optimize ("O3")

///////////////
// cvt calls might only work for values of magnitude less than 8388608?

// http://msdn.microsoft.com/en-us/library/yxty7t75.aspx
// http://sepwww.stanford.edu/data/media/public/sep//claudio/Research/Prst_ExpRefl/ShtPSPI/intel/cce/10.1.015/include/emmintrin.h
// http://sepwww.stanford.edu/data/media/public/sep//claudio/Research/Prst_ExpRefl/ShtPSPI/intel/cce/10.1.015/include/emm_func.h
// http://stackoverflow.com/questions/16447265/my-sse2-flooring-function-has-some-problems
// http://www.ibiblio.org/gferg/ldp/GCC-Inline-Assembly-HOWTO.html#s3
inline float math_floor(float a)
{
	__m128 a128 = _mm_set1_ps( a );
    __m128 trunc = _mm_cvtepi32_ps( _mm_cvttps_epi32( a128 ) );
    __m128 cmp = _mm_cmpgt_ps( trunc, a128 ); // check for -ve a
    __m128 and = _mm_and_ps( cmp, _mm_set1_ps( 1.0f ) );
    __m128 result = _mm_sub_ps( trunc, and );
    return _mm_cvtss_f32( result );
}

inline float math_ceil(float a)
{
	__m128 a128 = _mm_set1_ps( a );
    __m128 trunc = _mm_cvtepi32_ps( _mm_cvttps_epi32( a128 ) );
    __m128 cmp = _mm_cmplt_ps( trunc, a128 ); // check for +ve a
    __m128 and = _mm_and_ps( cmp, _mm_set1_ps( 1.0f ) );
    __m128 result = _mm_add_ps( trunc, and );
    return _mm_cvtss_f32( result );
}

inline float math_trunc(float a)
{
	__m128 a128 = _mm_set1_ps( a );
    __m128 trunc = _mm_cvtepi32_ps( _mm_cvttps_epi32( a128 ) );
    return _mm_cvtss_f32( trunc );
}

inline float math_frac(float a)
{
	__m128 a128 = _mm_set1_ps( a );
    __m128 trunc = _mm_cvtepi32_ps( _mm_cvttps_epi32( a128 ) );
    __m128 frac = _mm_sub_ps( a128, trunc );
    return _mm_cvtss_f32( frac );
}

inline float math_fmod(float a, float b)
{
    __m128 c = _mm_div_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) );
    __m128 c_trunc = _mm_cvtepi32_ps( _mm_cvttps_epi32( c ) );
    __m128 base = _mm_mul_ps( c_trunc, _mm_set1_ps( b ) );
    __m128 result = _mm_sub_ps( _mm_set1_ps( a ), base );
    return _mm_cvtss_f32( result );
}

static __attribute__((aligned (16))) unsigned long MATH_FABS_MASK = 0x7FFFFFFF;

#if 0 // math_fabs

inline float math_fabs(float a)
{
#if 0
	// 3.5 instr, 1 branch
	return a < .0f ? -a : a;
#else
	// 2 instr
	__m128 mask128 = _mm_load_ps( (float*)&MATH_FABS_MASK );
	__m128 pos128 = _mm_and_ps( _mm_set1_ps( a ), mask128 );
    return _mm_cvtss_f32( pos128 );
#endif
}

#else // math_fabs

#define math_fabs(a) \
	_mm_cvtss_f32( \
		_mm_and_ps( \
			_mm_set1_ps( (a) ), \
			_mm_load_ps( (float*)&MATH_FABS_MASK ) \
		) \
	)

#endif // math_fabs

#if 0 // math_sign

inline float math_sign(float a)
{
	return a < .0f ? -1.f : +1.f;
}

#else // math_sign

#define math_sign(a) \
	( (a) < .0f ? -1.f : +1.f )

#endif // math_sign

#if 0 // math_interp_wrap

// to test...
inline void math_interp_wrap(float* pa, float slope, float min, float max)
{
	// 26 instr, 1 branch
	float b = *pa + slope;
	float target = slope > 0.f ? max : min;
	float other = slope > 0.f ? min : max;
	float test = b - target;
	*pa = test *  math_sign( slope ) > 0.f ? other + test : b;
}

#else // math_interp_wrap

#define math_interp_wrap(pa, slope, min, max) \
	(*pa) = ( \
		(*pa) + (slope) < (min) \
		? (*pa) + (slope) + (max) - (min) \
		: ( \
			(*pa) + (slope) > (max) \
			? (*pa) + (slope) + (min) - (max) \
			: (*pa) + (slope) \
		) \
	)

#endif // math_interp_wrap

#if 0 // math_interp_pong

inline void math_interp_pong(float* pa, float* pslopesign, float slope, float min, float max)
{
	// 32 instr, 1 branch
	(*pa) += (*pslopesign) * (slope);
	(*pslopesign) *= math_fabs( *pa ) > (max) ? -1.f : +1.f;
 	if( math_fabs( *pa ) > (max) ) (*pa) += math_sign( *pa ) * ( (max) - math_fabs(*pa) );
}

#else // math_interp_pong

#define math_interp_pong(pa, pslopesign, slope, min, max) \
	do { \
		(*pa) += (*pslopesign) * (slope); \
		(*pslopesign) *= math_fabs( *pa ) > (max) ? -1.f : +1.f; \
		if( math_fabs( *pa ) > (max) ) (*pa) += math_sign( *pa ) * ( (max) - math_fabs(*pa) ); \
	} while(0)

#endif // math_interp_pong

#if 0 // math_interp_clamp

inline void math_interp_clamp(float* pa, float slope, float target)
{
#if 0
	// 5 instr, 1 branch
	float b = *pa + slope;
	*pa = slope > 0.f ? (b > target ? target : b) : (b < target ? target : b);
#else
	// 19 instr
	float b = *pa + slope;
	*pa = ( b - target ) *  math_sign( slope ) > 0.f ? target : b;
#endif
}

#else // math_interp_clamp

#define math_interp_clamp(pa, slope, target) \
	(*pa) = ( ( (*pa) + (slope) ) - (target) ) *  math_sign( (slope) ) > 0.f ? (target) : ( (*pa) + (slope) )

#endif // math_interp_clamp

// -1 to +1 over -pi/2 to pi/2 and repeating (looks like sin)
inline float math_saw(float a)
{
#if 0
	// 50 instr, 10 memy acc
	float frac = math_frac( a / M_PI_2 );
	float q1 = ( (Uint32)(a) + 1 ) % 2;
	float q2 = ( (Uint32)(a) ) % 2;
	return q1 * ( frac ) + q2 * ( frac - 1.f );
#else
	// 20 instr, 4 memy acc
	__m128 a128 = _mm_div_ps( _mm_set1_ps( a ), _mm_set1_ps( M_PI_2 ) );
	__m128i i128 = _mm_cvttps_epi32( a128 );
    __m128 trunc = _mm_cvtepi32_ps( i128 );
    __m128 frac = _mm_sub_ps( a128, trunc );
	__m128 quad2 = _mm_cvtepi32_ps( _mm_and_si128( i128, _mm_set1_epi32( 1 ) ) );
	__m128 quad1 = _mm_cvtepi32_ps( _mm_and_si128( _mm_add_epi32( i128, _mm_set1_epi32( 1 ) ), _mm_set1_epi32( 1 ) ) );
	__m128 saw2 = _mm_mul_ps( quad2, _mm_sub_ps( frac, _mm_set1_ps( 1.0f ) ) );
	__m128 saw1 = _mm_mul_ps( quad1, frac );
	__m128 result = _mm_add_ps( saw1, saw2 );
    return _mm_cvtss_f32( result );
#endif
}

inline float math_sin(float vlf_osc3_angle)
{
#if 0
	float result;
	asm( "fld %1; fsin; fstp %0;" : "=m" (result) : "m" (vlf_osc3_angle) );
	return result;
#else
	// only valid between -pi/2 and +pi/2
	// 19 instr, <30 cyc (core2), 4 reads
	const float B = 4.f / M_PI;
	const float C = -4.f / (M_PI * M_PI);
	const float P = 0.225f;

	float y = B * vlf_osc3_angle + C * vlf_osc3_angle * math_fabs( vlf_osc3_angle );
	return P * (y * math_fabs( y ) - y) + y;
#endif
}

#if 1 // math_percent_similar

// probably only safe for +ve integers...
inline float math_percent_similar(float a, float b)
{
#if 0
	// 8 instr, 3 br
	return (a == b) ? 1.f : a < b ? ( a / b ) : ( b / a );
#elif 0
	// 14 instr, 2 br
	return (a == b) ? 1.f : (a < b ? a : b) / (a > b ? a : b);
#else
	// 17 instr
	__m128 min = _mm_min_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) );
    __m128 min_zero = _mm_and_ps( _mm_cmpeq_ps( min, _mm_set1_ps( .0f ) ), _mm_set1_ps( 1.f ) );
	__m128 max = _mm_max_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) );
    __m128 max_zero = _mm_and_ps( _mm_cmpeq_ps( max, _mm_set1_ps( .0f ) ), _mm_set1_ps( 1.f ) );
	__m128 result = _mm_div_ps( _mm_add_ps( min, min_zero ), _mm_add_ps( max, max_zero ) );
    return _mm_cvtss_f32( result );
#endif
}

#else // math_percent_similar

#define math_percent_similar(a,b) \
    _mm_cvtss_f32( \
		_mm_div_ps( \
			_mm_add_ps( \
				_mm_min_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) ), \
				_mm_and_ps( _mm_cmpeq_ps( _mm_min_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) ), _mm_set1_ps( .0f ) ), _mm_set1_ps( 1.f ) ) \
			), \
			_mm_add_ps( \
				_mm_max_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) ), \
				_mm_and_ps( _mm_cmpeq_ps( _mm_max_ps( _mm_set1_ps( a ), _mm_set1_ps( b ) ), _mm_set1_ps( .0f ) ), _mm_set1_ps( 1.f ) ) \
			) \
		) \
	)

#endif // math_percent_similar

static unsigned long math_rnd_seed = 0;
inline void math_rand_seed()
{
	math_rnd_seed = perf_tsc();
}
inline unsigned short math_rand_u16()
{
	math_rnd_seed = math_rnd_seed * 1566083941 + 69069; // per Knuth, Table 1 section 3.3.4 TAOCP vol 2 1981
	return math_rnd_seed >> 16; // most entropy in high bits
}

#define MATH_RAND_U16() \
	( ( math_rnd_seed = math_rnd_seed * 1566083941 + 69069 ) >> 16 )

#if 1 // math_rand

float math_rand()
{
#if 0
	// 8 instr
	return (float)math_rand_u16() * 0.00001525902f;
#else
	// 4 instr
	__m128 inv65535 = _mm_set1_ps( 0.00001525902f /* 1 / 65535 */ );
	__m128 randShort = _mm_cvtepi32_ps( _mm_set_epi16( 0, 0, 0, 0, 0, 0, 0, math_rand_u16() ) );
	__m128 result = _mm_mul_ps( randShort, inv65535 );
    return _mm_cvtss_f32( result );
#endif
}

#else // math_rand

#define math_rand() \
	_mm_cvtss_f32( \
		_mm_mul_ps( \
			_mm_cvtepi32_ps( _mm_set_epi16( 0, 0, 0, 0, 0, 0, 0, MATH_RAND_U16() ) ), \
			_mm_set1_ps( 0.00001525902f /* 1 / 65535 */ ) \
		) \
	)

#endif // math_rand

#pragma GCC pop_options

///////////////////////

int main(int argc, char **argv)
{
	math_rand_seed();

	if( SDL_Init( SDL_INIT_AUDIO | SDL_INIT_VIDEO ) != 0 )
	{
		printf("SDL_Init failed\n");
		exit(EXIT_FAILURE);
	}

	if( SDL_SetVideoMode( 320, 160, 0, SDL_HWSURFACE ) == NULL )
	{
		printf("SDL_SetVideoMode failed\n");
		exit(EXIT_FAILURE);
	}
	SDL_WM_SetCaption( "A Simple SDL Window", 0 );

	if( SDL_OpenAudio(&sdl_configrequest, &sdl_config) != 0 )
	{
		printf("SDL_OpenAudio failed\n");
		exit(EXIT_FAILURE);
	}
	printf("SDL_OpenAudio samples %d\n", sdl_config.samples);
	printf("SDL_OpenAudio size %d\n", sdl_config.size);

	buf_init();

	synth_renderEvent( -1, 0.f, 0.f, 0.f );

	/////////////////

	wii_wiimotes =  wiiuse_init(MAX_WIIMOTES);
	int found = wiiuse_find(wii_wiimotes, MAX_WIIMOTES, 2 /* seconds to wait */ );
	if( !found )
	{
		printf("No wiimotes found.\n");
		exit(EXIT_FAILURE);
	}

	wii_connected = wiiuse_connect(wii_wiimotes, MAX_WIIMOTES);
	if( !wii_connected )
	{
		printf("Failed to connect to any wiimote.\n");
		exit(EXIT_FAILURE);
	}
	printf("Connected to %i wiimotes (of %i found).\n", wii_connected, found);

	wiiuse_motion_sensing(*wii_wiimotes, 1);

	///////////////

	printf("Creating SDL thread...");
    synth_threadID = SDL_CreateThread(synth_renderThread, (void *)NULL);
    if( NULL == synth_threadID )
    {
        printf("\nSDL_CreateThread failed: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
    }
	printf("made\n");

	printf("Waiting for thread to start...");
	BUF_WAIT_FIRSTUSE();
	printf("started\n");

	SDL_PauseAudio(0);

	SDL_Event event;
	while( !app_stop )
	{

		if( wiiuse_poll(wii_wiimotes, MAX_WIIMOTES) && wii_wiimotes[0]->event == WIIUSE_EVENT )
		{
			float roll = wii_wiimotes[0]->orient.roll, pitch = wii_wiimotes[0]->orient.pitch, yaw = wii_wiimotes[0]->orient.yaw;

			SPINLOCK_WAIT( synth_renderSpinlock );

			Uint8 u8event = -2; // only accel event
			if( IS_JUST_PRESSED( wii_wiimotes[0], WIIMOTE_BUTTON_A ) ) u8event = -3;
			else if( IS_HELD( wii_wiimotes[0], WIIMOTE_BUTTON_A ) ) u8event = -4;
			else if( IS_RELEASED( wii_wiimotes[0], WIIMOTE_BUTTON_A ) ) u8event = -5;

			synth_renderEvent( u8event, roll, pitch, yaw );

#if 0

			static Uint32 wmtc = 0-1;
			if( wmtc == 0-1 ) { wmtc = SDL_GetTicks() + 500; }
			if( SDL_GetTicks() > wmtc )
			{
				wmtc = SDL_GetTicks() + 500;
				printf("wiimote R %+03.1f P %+03.1f Y %+03.1f\n", roll, pitch, yaw);
			}

#endif

		}

		int evt = 0;
		while( SDL_PollEvent( &event ) )
		{
			evt++;
			app_stop = (event.type == SDL_QUIT) | ((event.type == SDL_KEYDOWN) & (event.key.keysym.sym == SDLK_ESCAPE));
			if( event.type == SDL_KEYDOWN )
			{
				SPINLOCK_WAIT( synth_renderSpinlock );
				synth_renderEvent( event.key.keysym.sym, 0.f, 0.f, 0.f );

				printf("key %x\n", event.key.keysym.sym);
			}
		}

#if defined(ENABLE_RENDERTIMING)

		static Uint32 tick = 0;
		if( (++tick) % 800 == 0 )
		{
			double avgRenderCycles = (double)synth_loopCycles / (double)synth_loopIter;
			printf("avg render cycles %f, avg render speed %f hz\n", avgRenderCycles, (double)2.16e+9 / avgRenderCycles );
		}

#endif // ENABLE_RENDERTIMING

	}

	printf("Waiting for thread to finish...");
	SDL_WaitThread(synth_threadID, &synth_threadCode);
	printf("finished\n");

    SDL_CloseAudio();
	SDL_Quit();

	buf_cleanup();

	wiiuse_cleanup(wii_wiimotes, MAX_WIIMOTES);

	return 0;
}

///////////////////////

// the optimizations don't seem to work with a multithreaded sdl app, so we turn them on locally only
#pragma GCC push_options
#pragma GCC optimize ("O3")

FLT_INTERP_NEW( twist_reference );
FLT_INTERP_NEW( twist_follow );
FLT_FUNC_NEW( twist_similar );
FLT_FUNC_NEW( twist_filtered );

FLT_SWEEP_NEW( vlf_osc1_angle );
FLT_FUNC_NEW( vlf_osc1 );

FLT_SWEEP_NEW( vlf_osc2_angle );
FLT_FUNC_NEW( vlf_osc2 );

FLT_INTERP_NEW( vlf_osc3_angleTarget );
FLT_SWEEP_NEW( vlf_osc3_angle );
FLT_FUNC_NEW( vlf_osc3 );

FLT_NEW( vlf_mix );

FLT_INTERP_NEW( roll_reference );
FLT_INTERP_NEW( roll_follow );
FLT_FUNC_NEW( roll_similar );
FLT_FUNC_NEW( roll_filtered );

FLT_SWEEP_NEW( mf_osc1_angle );
FLT_FUNC_NEW( mf_osc1 );

FLT_SWEEP_NEW( mf_osc2_angle );
FLT_FUNC_NEW( mf_osc2 );

FLT_IIRF_NEW( spark_iirf );

FLT_SWEEP_NEW( spark_osc1_angle );
FLT_FUNC_NEW( spark_osc1 );

FLT_INTERP_NEW( spark_osc2_ampTarget );
FLT_SWEEP_NEW( spark_osc2_angle );
FLT_FUNC_NEW( spark_osc2 );

FLT_FUNC_NEW( spark_noise );

FLT_NEW( spark_mix );

FLT_NEW( mix );

float vlf = 150.f;
float mf = 280.f;
float hf = 440.f;

void synth_renderEvent( Sint8 code, float rollIn, float pitchIn, float yawIn )
{
	float k1HzDelta = 2.f * M_PI / (float)sdl_config.freq;
	float twist = 99.f * ( math_fabs( pitchIn ) + math_fabs( rollIn ) ) / 360.f + 1.f /* avoid 0 for similar calculation */ ; // [1,100]
	float roll = 99.f * math_fabs( rollIn ) / 180.f + 1.f /* avoid 0 for similar calculation */ ; // [1,100]

	if( code == -1 ) // init
	{
		FLT_SET( twist_reference, 1.f );
		FLT_SET( twist_follow, 1.f );

		FLT_SET( roll_reference, 1.f );
		FLT_SET( roll_follow, 1.f );

		FLT_SWEEP_PONG_CONFIG( vlf_osc1_angle, vlf * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SWEEP_PONG_CONFIG( vlf_osc2_angle, .25f * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SWEEP_PONG_CONFIG( vlf_osc3_angle, vlf * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SET( vlf_osc3_angleTarget, vlf * k1HzDelta );

		FLT_SWEEP_PONG_CONFIG( mf_osc1_angle, 8.f * mf * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SWEEP_PONG_CONFIG( mf_osc2_angle, mf * k1HzDelta, -M_PI_2, +M_PI_2 );

		FLT_SET( spark_iirf, 0.f );
		FLT_SWEEP_PONG_CONFIG( spark_osc1_angle, .25f * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SWEEP_PONG_CONFIG( spark_osc2_angle, hf * k1HzDelta, -M_PI_2, +M_PI_2 );
		FLT_SET( spark_osc2_ampTarget, 0.f );
	}
	else
	{
		FLT_INTERP_TRIGGER( twist_reference, 10.f * k1HzDelta /* fast to avoid click */, twist );
		FLT_INTERP_TRIGGER( twist_follow, 1.f * k1HzDelta, twist );

		FLT_INTERP_TRIGGER( roll_reference, 10.f * k1HzDelta /* fast to avoid click */, roll );
		FLT_INTERP_TRIGGER( roll_follow, 1.f * k1HzDelta, roll );

		FLT_INTERP_TRIGGER( vlf_osc3_angleTarget, .05f * k1HzDelta, vlf * k1HzDelta );

		switch( code )
		{
			case -3 /* press */: FLT_INTERP_TRIGGER( spark_osc2_ampTarget, 10.f * k1HzDelta /* fast to avoid click */, 1.f ); break;
			case -4 /* hold */: if( math_rand() < .1f ) { FLT_IIRRF_TRIGGER( spark_iirf, (float)sdl_config.freq, 1000.f, .0001f, .0001f ); } break;
			case -5 /* release */: FLT_INTERP_TRIGGER( spark_osc2_ampTarget, 10.f * k1HzDelta /* fast to avoid click */, 0.f ); break;
		}
	}
}

float synth_render( Uint32 uT  )
{
	float k1HzDelta = 2.f * M_PI / (float)sdl_config.freq;

	FLT_INTERP( twist_reference );
	FLT_INTERP( twist_follow );
	FLT_FUNC_SIMILAR( twist_similar, twist_reference, twist_follow );
	FLT_MIX2_CROSS( twist_filtered, .90f, twist_filtered, twist_similar );

	FLT_INTERP( roll_reference );
	FLT_INTERP( roll_follow );
	FLT_FUNC_SIMILAR( roll_similar, roll_reference, roll_follow );
	FLT_MIX2_CROSS( roll_filtered, .90f, roll_filtered, roll_similar );

	////////

	FLT_SWEEP_PONG( vlf_osc1_angle );
	FLT_FUNC_SIN( vlf_osc1, 1.f, FLT_OUTPUT( vlf_osc1_angle ) );

	FLT_SWEEP_PONG( vlf_osc2_angle );
	FLT_FUNC_SIN( vlf_osc2, 1.f, FLT_OUTPUT( vlf_osc2_angle ) );
	
	FLT_INTERP( vlf_osc3_angleTarget );
	FLT_SWEEP_PONG_RATED( vlf_osc3_angle, FLT_OUTPUT( vlf_osc3_angleTarget ) * ( 1.f + .1f * FLT_OUTPUT( twist_filtered ) ) );	
	FLT_FUNC_SIN( vlf_osc3, .8f + .2f * FLT_OUTPUT( vlf_osc2 ), FLT_OUTPUT( vlf_osc3_angle ) );
	
	FLT_MIX2( vlf_mix, vlf_osc3, vlf_osc1 );

	////////
	
	FLT_SWEEP_PONG( mf_osc1_angle );
	FLT_FUNC_SAW( mf_osc1, 1.f, FLT_OUTPUT( mf_osc1_angle ) );
	
	FLT_SWEEP_PONG( mf_osc2_angle );
	FLT_FUNC_SIN( mf_osc2, ( 1.f - FLT_OUTPUT( roll_filtered ) ) * ( .5f + .1f * FLT_OUTPUT( mf_osc1 ) ), FLT_OUTPUT( mf_osc2_angle ) );

	////////

	FLT_SWEEP_PONG( spark_osc1_angle );
	FLT_FUNC_SIN( spark_osc1, 1.f, FLT_OUTPUT( spark_osc1_angle ) );

	FLT_INTERP( spark_osc2_ampTarget );
	FLT_SWEEP_PONG_RATED( spark_osc2_angle, hf * k1HzDelta * ( 1.f + .1f * FLT_OUTPUT( twist_filtered ) ) );
	FLT_FUNC_SIN( spark_osc2, ( .65f + .25f * FLT_OUTPUT( spark_osc1 ) ) * FLT_OUTPUT( spark_osc2_ampTarget ), FLT_OUTPUT( spark_osc2_angle ) );

	FLT_IIRRF( spark_iirf );
	FLT_FUNC_WHITE( spark_noise, FLT_OUTPUT( spark_iirf ) );

	FLT_MIX2( spark_mix, spark_osc2, spark_noise );

	////////

	FLT_MIX3( mix, vlf_mix, mf_osc2, spark_mix );
	return FLT_OUTPUT( mix );
}

#pragma GCC pop_options

