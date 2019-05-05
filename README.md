# WARNING

This library is under active development.  Do not hitch a wagon to it because it can and will change drastically.

# RTSeis

RTSeis is a toolkit intended to aid in real-time filtering of seismic signals.  As real-time signals processing is a superset of post-processing, RTSeis naturally has a feature-rich suite of post-processing algorithms.  While the core library is written in C++, note that there exists Python bindings built around the core library that expose a subset of the library's functionality.  These Python wrappers can outperform NumPy and SciPy equivalent implementations by a factor of 2-4 and outperform ObsPy by factors of 10 to 100.

To get started go to the directory of your choosing and get the code

    git clone https://github.com/bakerb845/rtseis.git
    
Next, descend into rtseis and build the documentation

    cd rtseis
    doxygen Doxyfile

Finally, open the documentation with a web browser.  This should provide build examples and usage instructions.  For GCC users, an always valid build example is in .circleci/config.yml.
