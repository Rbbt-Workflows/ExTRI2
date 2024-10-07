import json
import re
from IPython.display import display, HTML
from IPython.display import Markdown


# Create markdown shortcuts
def md(text):
    display(Markdown(text))

def display_header(level):
    '''md(f'<hX>{text}</hX}>') -> hX(text) '''
    def header_function(text):
        md(f'<h{level}>{text}</h{level}>')
    return header_function

h1 = display_header(1)
h2 = display_header(2)
h3 = display_header(3)
h4 = display_header(4)
h5 = display_header(6)
h6 = display_header(6)

def bold(text):
    md(f'<b>{text}</b>')

## OTHER USEFUL FUNCTIONS FOR THE NOTEBOOK

def highlight_words(text: str, words: list, wrapper: str = "strong"):
    '''Replace occurrences of a list of words in a text with the words wrapped in <strong> tags'''
    for word in words:
        text = text.replace(word, f'<{wrapper}>{word}</{wrapper}>')
    text = text.replace('.', '.<br>') 
    return text

def table_from_dict(title, d: dict, heading='h3'):
    '''Create a HTML table from a dictionary'''
    html_str = f'<{heading}>{title}</{heading}><table>'
    for key, value in d.items():
        html_str += f'<tr><td><b>{key}</b></td><td>{value}</td></tr>'
    html_str += '</table>'
    display(HTML(html_str))


def table_of_contents(notebook_path):
    '''Displays the Table of contents of a given notebook path'''
    toc = []
    
    with open(notebook_path, "r") as f:
        cells = json.load(f)["cells"]
    
    for cell in cells:
        if cell["cell_type"] == "markdown":
            for line in cell["source"]:
                match = re.search("^#+ \w+", line)
                if match:
                    level = len(line) - len(line.lstrip("#"))
                    link = line.strip(" #\n").replace(" ", "-")
                    toc.append(
                        2 * (level - 2) * " "
                        + (("- [") if level > 1 else '\n[')
                        + line.strip(" #\n")
                        + f"](#{link})"
                    )
    h3('Table of contents')
    md('\n'.join(toc))