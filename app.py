import dash
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import app_set as setup
import numpy as np
import time 
import tkinter as tk
root = tk.Tk()
height_px = root.winfo_screenheight()

t1 = time.time()

# preporcess data
data = setup.Preprocess(hdf5_output_config_file='output_GJ.hdf5',
                        hdf5_input_config_file='network-synapses.hdf5')
t2 = time.time()
print(t2 - t1)
print('data loaded!')


# parallel coordinate graph--------------------------------------------------

fig = go.Figure(data=
    go.Parcoords(
        line = dict(color = data.df['has_db'],
                   colorscale = [[0,'green'],[1,'red']]),
        dimensions = data.dim, name='pkey',
        labelfont={'size':20}
    )
)


# gap junctions ------------------------------------------------------------

# compare #gj and mean somatic distance etc of cells with and without DB
fig3 = ff.create_distplot([ data.gj_count_for_hist['neg'],
                            data.gj_count_for_hist['pos']], 
                            ['neg','pos'], 
                            bin_size=4,
                            colors=['green','red'],
                            rug_text=[  data.gj_count_for_hist['rug_text']['neg'],
                                        data.gj_count_for_hist['rug_text']['pos']])
fig3['layout'].update(height=int(height_px/2160*350))



app = dash.Dash(__name__)

# layout: html/css setup ---------------------------------------------------------------

app.layout = html.Div(children=[
    # first row: title
    html.H1(children='Distribution of optimization parameters'),
    
    # second row: parallel coordinate plot
    html.Div([
        dcc.Graph(figure=fig)
        ], style={'width': '90%', 'display': 'inline-block', 'padding': '0 20'}),
    
    # number of gapjunctions in cells with vs without depolarization block
    html.Div([
        dcc.Graph(figure=fig3) 
        ], style={'width': '30%', 'display': 'inline-block', 'padding': '0 20'}),
])    
#    # third row: voltage traces of depol block
#    html.Div([
#        dcc.Graph(figure=fig2)
#        ], style={'width': '90%', 'display': 'inline-block', 'padding': '0 20'}),
    


# ---------------------------------------------------------------------------------------
'''
@app.callback(
    dash.dependencies.Output('cells_with_DB', 'figure'),
    [dash.dependencies.Input('cell_selecter', 'value')])
def update_voltage_DB(cellid):
    
    # all patterns
    patterns    = mapper[data_set][pattern_size]['patterns' ]
    
    df          = DF[data_set]
    filtered_df = df[df.pattern_len == pattern_size]
    traces      = []
    
    for i in filtered_df.Model.unique():
        x       = []
        y       = []
        text    = []
        for j in filtered_df.pattern.unique():
            df_by_model = filtered_df[  (filtered_df['Model']==i)   &
                                        (filtered_df['DA']==da)     &
                                        (filtered_df['pattern']==j) ]
            # get y-value (last in list)
            y.append( df_by_model['count'].tolist()[-1] )
            
            # get section (last in pattern) and somatic distance
            s     = str(patterns[j][-1])  
            x.append( morphology['stat'][int(s)][x_axis_key] )
            text += ['sec:'+s+' pattern:'+str(j)+' stem:'+str(sec2stem[s])]
          
        traces.append(go.Scatter(
            x=x,
            y=y,
            text=text,
            mode='markers',
            opacity=0.5,
            showlegend=True,
            marker={    'color': colors[i%len(colors)],
                        'size': 10,
                        'line': {'width': 1, 'color': 'white'}
            },
            name=str(i)
        ))
    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'type': 'linear', 'title': x_axis_key},
            yaxis={'title': 'count'},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 1, 'y': 0},
            hovermode='closest',
            height=350
    )}
'''
      
# ---------------------------------------------------------------------------------------
'''
# somatic Vm (of cells with DB) -----------------------------------------------
traces      = []     
for cellid,item in data.cells_with_db.items():
    # full trace
    traces.append(go.Scatter(
        x=data.time,
        y=item['v'],
        text=cellid,
        mode='lines',
        opacity=0.5,
        showlegend=True,
        legendgroup=cellid,
        marker={ 'color': 'grey'},
        name=str(cellid)
    ))
    
    # depolarization blocks (DB)
    x = []
    y = []
    for i,val in enumerate(item['start_indx']):
        e = data.time[item['end_indx'][i]]
        m = np.max(item['v'][val:item['end_indx'][i]])
        x += [data.time[val],e]
        y += [m,m]
        x.append(None)
        y.append(None)
    
    traces.append(go.Scatter(
        x=x,
        y=y,
        text='db:{}'.format(cellid),
        mode='lines',
        opacity=0.5,
        showlegend=False,
        legendgroup=cellid,
        marker={    'color': 'red' },
        name='{}:db'.format(cellid)
    ))    
            
fig2 = {
        'data': traces,
        'layout': go.Layout(
            xaxis={'type': 'linear', 'title': 'time'},
            yaxis={'title': 'Vm'},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 1, 'y': 0},
            hovermode='closest',
            height=350
    )}
'''
    
if __name__ == '__main__':
    app.run_server(debug=True)
