# outline for a myst_nb project with sphinx
# build with: sphinx-build -nW --keep-going -b html . ./_build/html

# Make sure autodoc can find code
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

# load extensions
extensions = ["myst_nb", "sphinx_rtd_theme", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]


autodoc_class_signature = "separated"
autodoc_member_order = "bysource"

myst_enable_extensions = ["dollarmath"]

html_logo = "_static/logo.png"

# specify project details
master_doc = "index"
project = "WAM2layers"

# basic build settings
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
nitpicky = True

html_theme = "sphinx_rtd_theme"
