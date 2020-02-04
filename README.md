# RTSeis

RTSeis is a toolkit intended to aid in real-time filtering of seismic signals.  As real-time signals processing is a superset of post-processing, RTSeis naturally has an extensive suite of post-processing algorithms.  While the core library is written in C++, it is worth nothign that there exist Python bindings that expose a subset of the library's functionality.  These Python wrappers can outperform NumPy/SciPy equivalent implementations by 2-4x and ObsPy by factors in excess of 10x.

To get started go to the directory of your choosing and get the code

    git clone https://github.com/uofuseismo/rtseis.git
    
Next, descend into rtseis and build the documentation

    cd rtseis
    doxygen Doxyfile

Finally, open the documentation with a web browser.  This should provide build examples and usage instructions.  For GCC users, an always valid build example is in .circleci/config.yml.

## Disclaimers

The API is still being refined but the usage is basically worked.  For the low-level stuff you would do something

   1.  Initialize the filter
   2.  If needed, set some initial conditions.
   3.  Continually apply the data.
   4.  For real-time processing reset the filter to its initial conditions in the presence of a gap.
   5.  Clear the filter when you are finished (the destructor will automatically release all memory as well).

More straightforward are higher level modules like:

    postProcessing/singleChannel/waveform.hpp

The only potential major API change I foresee would be making the low-level real-time and post-processing behavior a filter specific to a template class.  This means that instead of doing

    filter.initialize(processing_mode);

you would instead do something like, 

    Filter<precision, processing\_mode> filter;

Again, if you use the high-level API you won't know this switch was made.  Therefore, I'd recommend sticking with the high-level API since it is considerably less likely to change in any significant way.

# Build Status:

[![CircleCI](https://circleci.com/gh/uofuseismo/rtseis.svg?style=svg&circle-token=build_status)](https://circleci.com/gh/uofuseismo/rtseis)
 
