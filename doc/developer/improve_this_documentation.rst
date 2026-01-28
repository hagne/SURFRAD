Help improve this documentation
================================

Build this documentation locally 
--------------------------------
Requirements
^^^^^^^^^^^^
To build the documentation locally, you need to have Python installed along with the Sphinx documentation generator
and the necessary dependencies. The Read the Docs theme is not bundled with Sphinx, so install ``sphinx-rtd-theme``.
You can install Sphinx and the required packages using conda or pip.
Using conda:
   conda install sphinx sphinx-rtd-theme
Using pip:
   pip install sphinx sphinx-rtd-theme

Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^
Once you have the requirements installed, navigate to the root directory of the surfradpy documentation (e.g. doc/) and run:
    make html
This will generate the HTML documentation in the _build/html directory.
You can then open the index.html file in a web browser to view the documentation.
