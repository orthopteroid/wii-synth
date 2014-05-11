wii-synth
=========

Linux-platform gesture-based synthesiser designed around using a wii controller as an input device. Uses SDL and WiiUse libraries and a C-based SSE2 assembly-framework to build the synth using GCC. This synth is really a toy, as the synth network is placed in a static memory using hand-coded using C macros. But it does show how to do low-latency, event-driven realtime audio rendering. It uses some x86 and SSE assembler (actually lots of that really).

## Oscillators

The heart of the synth is how it's oscillators work. There are two types: the infinite-impulse-response-function (an ODE solution to the spring/damper oscillator, as decribed inPerry Cook's "Real Sound Synthesis for Interactive Applications", p43) and a linear-angle-sweep method (LAS). The latter is mainly used for the effects in the synth network so I'll describe that here.

The LAS method simply adds an increment to an angle value and restricts the result to the interval between -pi/2 and +pi/2. In order to restrict the value the LAS value is either 'ponged' (bounces back from one of the extremes back into the interval) or 'wraped' (the amount the value exceeds past one extreme is used as an offest to the other extreme into the interval). Restricting the LAS between -pi/2 and +pi/2 has the advantage of facilitating fast & accurate approximations to trig functions that work well just on that interval (see http://devmaster.net/posts/9648/fast-and-accurate-sine-cosine). Additionally, other functions can be written that only need to function on that interval: like functions for triangle and sawtooth waveforms that can change timbre by introducing other frequencies.

## Synth Network

The synth network code is made of three parts: the declaration, the event handler and the render loop. This doesn't count the framework that's needed to use SDL for audio or the wiiuse library to watch for changes to the wii controller (those parts make up about 700 lines, or 70%, of the code here). All code is C99, using some simple gcc extensions (like `__COUNTER__`).

The synth network is built using macros that place the network elements in static memory at file-scope. Everything goes into `float synth_filterNetwork [];` Some network elements have state or parameters but all elements return a value - it works kind of like a functional programming network. For each element declared in the network it may claim several 'registers' from this static array. Managing all that static storage is done with the `FLT_???_NEW` macros which use `__COUNTER__` to do that. Offset 0 for any network element is accumed to the the output of that element, but additional outputs can be present (hmmm, I think that was an old iteration of the code). In any case, the network is stored in the static array.

### Events

The network's event handler (`synth_renderEvent`) initializes and makes changes to the calculations in the network in a way that is safe to do without interfering with the smooth construction of waveforms. For example, if a user event is suppose to change the pitch of a voice that is being rendered the pitch change should happen without any audible click as the waveform is adjusted. This kind of thing is accomplished using a network element called an interpolator that is programmed to adjust it's output value at a specific rate until it's output tracks the specified input. When an event comes in the handler tells the interpolator to target towards a different value. The actual interpolation is performed in the next step: the renderer.

The event handler takes an event 'code' that indicates what kind of event had happened: -1 is initialize and -3 to -5 are press, hold and release of the wii button. These values are mainly historial, as the code is a leftover from the polysynth project where it previously indicates a keyboard keypress.

### Rendering

The actual meat of the synth (`synth_render`) performs waveform rendering as well as making adjustments to the network interpolators that facilitate the changing of pitch and amplification from the user events. This function and what it calls to render sound is designed without the need of any branching (ie no if or decision statements). It, the `FLT_???` macros and the `math_???` functions they call have been specifically been designed together this way in order to render as fast as possible (ie realtime) by preventing cpu pipeline stalls. Of course, this code was developed on old Core2-Duo T7400 laptop so YMMV.

## Higher Level

Looking at the code from a higher level, it uses wiiuse and SDL, the latter of which starts a playback thread which plays back from a 4k buffer. The synth renderer runs in another thread and uses spinlocks and double-buffering in order to keep the SDL's playback thread fed and happy. spinlocks are also used to ensure that when events are handled (from the program main thread) the renderer is between buffers. You can see much of that logic in `synth_renderThread` but also look for the variable `synth_renderSpinlock ` which shows how waiting is done when a wii event occurs.

Cheers.
