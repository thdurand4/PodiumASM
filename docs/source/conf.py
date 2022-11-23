import podiumASM
from podiumASM.snakeWrapper.global_variable import *

# The short X.Y version.
version = podiumASM.__version__
# The full version, including alpha/beta/rc tags
release = podiumASM.__version__

rst_prolog = f"""
.. |tools_path| replace:: {GIT_TOOLS_PATH}
"""





# -- Project information -----------------------------------------------------
# General information about the project.
project = 'PodiumASM'
copyright = '2022, T Durand (CIRAD), S Bache (CIRAD), S Ravel (CIRAD)'
github_doc_root = 'https://github.com/thdurand4/PodiumASM/tree/master/docs/'
issues_github_path = 'https://github.com/thdurand4/PodiumASM/issues'

latex_authors = '''
Theo Durand (CIRAD),\\\\
Simon Bache (CIRAD),\\\\
Sebastien Ravel (CIRAD)\\\\
'''

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
              'sphinx.ext.viewcode',
              'sphinx.ext.autosectionlabel',
              'sphinx_copybutton',
              'sphinx_rtd_theme',
              'sphinx_click'
              ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = ['.rst', "md"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

master_doc = 'index'

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'analytics_id': 'UA-172723859-1',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'cyan',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': False,
    'titles_only': False
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "PodiumASM"

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "PodiumASM"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_images/PodiumASM_logo.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_images/PodiumASM_logo2.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True
