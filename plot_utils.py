from io import BytesIO
import base64
import matplotlib.pyplot as plt
import numpy as np
from wordcloud import WordCloud

def abstracts_wordcloud(abstracts, ax=None):
    """ take a list of abstract texts and output a wordcloud """ 
    wordcloud = WordCloud(max_font_size=50, 
                          max_words=100, 
                          mask=np.zeros((500,500), dtype=int),
                          background_color="white").generate(' '.join(abstracts))
    if ax is None:
        f = plt.figure(figsize=(12,8));
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.axis("off")
        return f
    else:
        ax.imshow(wordcloud, interpolation="bilinear")


def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', bbox_inches='tight', pad_inches=0, **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)