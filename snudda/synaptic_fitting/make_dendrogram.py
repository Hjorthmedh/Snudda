from neurom.view.view import plot_dendrogram
from neurom.view.common import get_figure, plot_style
from neurom import load_neuron
from neurom.view.dendrogram import Dendrogram, get_size, layout_dendrogram, move_positions
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection, PatchCollection
import pathlib
import matplotlib.pyplot as plt

from neurom import NeuriteType
TREE_COLOR = {NeuriteType.basal_dendrite: 'red',
              NeuriteType.apical_dendrite: 'purple',
              NeuriteType.axon: 'blue',
              NeuriteType.soma: 'red'}

def _get_dendrogram_shapes(dendrogram, positions, show_diameters):
    """Generates drawable patches for dendrogram.
    Args:
        dendrogram (Dendrogram): dendrogram
        positions (dict of Dendrogram: np.array): positions xy coordinates of dendrograms
        show_diameter (bool): whether to draw shapes with diameter or as plain lines
    Returns:
        List of matplotlib.patches.
    """
    color = TREE_COLOR[dendrogram.neurite_type]
    start_point = positions[dendrogram]
    end_point = start_point + [0, dendrogram.height]
    if show_diameters:
        shapes = [_as_dendrogram_polygon(dendrogram.coords + start_point, color)]
    else:
        shapes = [_as_dendrogram_line(start_point, end_point, color)]
    for child in dendrogram.children:
        shapes.append(_as_dendrogram_line(end_point, positions[child], color))
        shapes += _get_dendrogram_shapes(child, positions, show_diameters)
    return shapes

def _get_dendrogram_legend(dendrogram):
    """Generates labels legend for dendrogram.
    Because dendrogram is rendered as patches, we need to manually label it.
    Args:
        dendrogram (Dendrogram): dendrogram
    Returns:
        List of legend handles.
    """
    def neurite_legend(neurite_type):
        return Line2D([0], [0], color=TREE_COLOR[neurite_type], lw=2, label=neurite_type.name)

    if dendrogram.neurite_type == NeuriteType.soma:
        handles = {d.neurite_type: neurite_legend(d.neurite_type)
                   for d in [dendrogram] + dendrogram.children}
        return handles.values()
    return [neurite_legend(dendrogram.neurite_type)]

def _as_dendrogram_polygon(coords, color):
    return Polygon(coords, color=color, fill=True)
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

def _as_dendrogram_line(start, end, color):
    return FancyArrowPatch(start, end, arrowstyle='-', color=color, lw=2, shrinkA=0, shrinkB=0)

def make_dendrogram(path_to_neuron,fig_spec={'format' : 'tiff'}):


    nrn = load_neuron(path_to_neuron)

    path_neuron = pathlib.Path(path_to_neuron)

    out_dir = pathlib.Path('dendrogram',*path_neuron.parts[1:-1])

    out_dir.mkdir(parents=True, exist_ok=True)

    out_name = '_'.join([path_neuron.name.split('.')[0],'dendrogram'])

    name = str(out_dir / out_name)
    
    plot_dend(nrn,name,fig_spec)
    

def plot_dend(nrn,name,fig_spec=None):
    
    dendrogram = Dendrogram(nrn)
    positions = layout_dendrogram(dendrogram, np.array([0, 0]))
    w, h = get_size(positions)
    positions = move_positions(positions, np.array([.5 * w, 0]))
    fig, ax = get_figure()
    ax.set_xlim([-.05 * w, 1.05 * w])
    ax.set_ylim([-.05 * h, 1.05 * h])
    ax.set_title('Morphology Dendrogram')
    ax.set_xlabel('micrometers (um)')
    ax.set_ylabel('micrometers (um)')
    shapes = _get_dendrogram_shapes(dendrogram, positions, show_diameters=False)
    ax.add_collection(PatchCollection(shapes, match_original=True))

    ax.set_aspect('auto')
    ax.legend(handles=_get_dendrogram_legend(dendrogram))
    plt.savefig('.'.join([name,fig_spec['format']]),dpi=300)
