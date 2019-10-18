# WARNING

This library is under active development.  The higher level API in postProcessing/singleChannel/waveform.hpp is stabilizing however the lower level utilities in utilities are still evolving.

# RTSeis

RTSeis is a toolkit intended to aid in real-time filtering of seismic signals.  As real-time signals processing is a superset of post-processing, RTSeis naturally has an extensive suite of post-processing algorithms.  While the core library is written in C++, it is worth nothign that there exist Python bindings that expose a subset of the library's functionality.  These Python wrappers can outperform NumPy/SciPy equivalent implementations by 2-4x and ObsPy by factors in excess of 10x.

To get started go to the directory of your choosing and get the code

    git clone https://github.com/uofuseismo/rtseis.git
    
Next, descend into rtseis and build the documentation

    cd rtseis
    doxygen Doxyfile

Finally, open the documentation with a web browser.  This should provide build examples and usage instructions.  For GCC users, an always valid build example is in .circleci/config.yml.
