# outline for a myst_nb project with sphinx
# build with: sphinx-build -nW --keep-going -b html . ./_build/html

# load extensions
extensions = ["myst_nb", "sphinx_rtd_theme", "sphinx.ext.autodoc", "sphinx.ext.napoleon", "autoapi.extension"]

# Options for autodoc
autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
import os, sys
current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../src"))
sys.path.insert(0, target_dir)
print(target_dir)

myst_enable_extensions = ["dollarmath"]

html_logo = "_static/logo.png"

# specify project details
master_doc = "index"
project = "WAM2layers"

# basic build settings
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
nitpicky = True

html_theme = "sphinx_rtd_theme"
