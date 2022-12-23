# outline for a myst_nb project with sphinx
# build with: sphinx-build -nW --keep-going -b html . ./_build/html

# load extensions
extensions = ["myst_nb", "sphinx_rtd_theme"]

myst_enable_extensions = ["dollarmath"]

html_logo = "_static/WAM_logo_v3.png"

# specify project details
master_doc = "index"
project = "WAM2layers"

# basic build settings
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
nitpicky = True

html_theme = "sphinx_rtd_theme"
