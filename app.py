import dash
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import app_set as setup
import numpy as np
import time, pickle
import tkinter as tk
root = tk.Tk()
height_px = root.winfo_screenheight()

t1 = time.time()

load_pickle = 1
if load_pickle:
    with open('app_data.pkl', 'rb') as h:
        data = pickle.load(h)
else:
    # preporcess data
    data = setup.Preprocess(hdf5_output_config_file='output_GJ.hdf5',
                            hdf5_input_config_file='network-synapses.hdf5')
    data.xx = None
    data.gap_junctions.xx = None
    with open('app_data.pkl', 'wb') as h:
        pickle.dump(data, h, protocol=pickle.HIGHEST_PROTOCOL)
    
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
#fig['layout'].update(height=int(height_px/2160*600))

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
    
    # second row: parallel coordinate plot --------------------
    html.Div([
        dcc.Graph(figure=fig)
        ], style={'width': '95%', 'padding': '0 20'}),
    
    # morphology and gap junctions ---------
    html.H1(children='Gap junctions and morphologies'),
    html.Div([
        html.Div([
            html.H3('cellid'),
            dcc.Input(
                id='cellid',
                placeholder='Enter a value...',
                value=1,
                type="number"
                )], style={'width':'50%', 'display': 'inline-block', 'padding': '10 20'}),
        html.Div([html.Hr()], style={'padding': '10 20'}),
        # morphology
        html.Div([
            dcc.Graph(id='graph-morph') 
            ], style={'width': '90%', 'padding': '10 20'})
        ], style={'width': '30%', 'display': 'inline-block', 'vertical-align':'top', 'padding': '10 20'}),
    
    # voltage trace...
    html.Div([
        html.Div([
            html.H3('Show voltage?'),
            dcc.RadioItems(
                    id='show-voltage',
                    options=[
                        {'label': 'On', 'value': 1},
                        {'label': 'Off','value': 0}
                        ],
                    value=0, 
                    style={'width':'50%', 'display': 'inline-block', 'padding': '10 20'})
                ]),
        html.Div([html.Hr()], style={'padding': '10 20'}),
        # voltage trace
        html.Div([
            dcc.Graph(id='graph-voltage') 
            ], style={'width': '90%', 'padding': '10 20'})
        ], style={'width': '60%', 'display': 'inline-block', 'vertical-align':'top', 'padding': '10 20'}),
    
    # number of gapjunctions in cells with vs without depolarization block
    html.Div([
        dcc.Graph(figure=fig3) 
        ], style={'width': '30%', 'display': 'inline-block', 'padding': '0 20', 'vertical-align':'bottom'})
])    
#    # third row: voltage traces of depol block
#    html.Div([
#        dcc.Graph(figure=fig2)
#        ], style={'width': '90%', 'display': 'inline-block', 'padding': '0 20'}),
    


# ---------------------------------------------------------------------------------------

@app.callback(
    dash.dependencies.Output('graph-morph', 'figure'),
    [dash.dependencies.Input('cellid', 'value')]
    )
def update_morph_graph(cellid):
    
    if not cellid:
        cellid=1
    if cellid < 0:
        cellid=0
        
    
    # from cellid: get family and mkey -> data
    fam = data.id2meta[cellid]['family']
    mkey = data.id2meta[cellid]['mkey']
    mdata = data.gap_junctions.morph_data[fam][mkey]
    
    traces              = []
    
    traces.append(go.Scatter3d(
        x = np.array(mdata['morph_with_none']['x']).astype(float),
        y = np.array(mdata['morph_with_none']['y']).astype(float),
        z = np.array(mdata['morph_with_none']['z']).astype(float),
        mode='lines',
        opacity=0.7,
        showlegend=False,
        line=dict( color='grey', width=3 )
        ))
    
    
    if cellid in data.gap_junctions.gj:
    
        # gj location
        x = []
        y = []
        z = []
        c = []
        w = []
        hovertext = []
        
        for key,l in data.gap_junctions.gj[cellid].items():
            
            if key == 0:
                x.append( 0 )
                y.append( 0 )
                z.append( 0 )
                c.append( 'blue' )
                w.append( len(l)*5.0 )
                hovertext.append('SOMA!, ngj={}'.format(len(l)))
                continue
            
            key = key-1
            
            # get coordinates
            if key not in mdata['sec_coordinates']:
                print(key, 'not in mdata...!!!')
                print(data.gap_junctions.gj[cellid])
                continue
            
            coordinates = mdata['sec_coordinates'][key]
            x.append( coordinates[0] )
            y.append( coordinates[1] )
            z.append( coordinates[2] )
            c.append( 'blue' )
            w.append( len(l)*5.0 )
            hovertext.append('sec={}, ngj={}'.format(key,len(l)))
        
        traces.append(go.Scatter3d(
            x=x,
            y=y,
            z=z,
            text=hovertext,
            mode='markers',
            showlegend=False,
            marker=dict( size=w, color=c )
            ))
    # soma
    traces.append(go.Scatter3d(
        x=[0],
        y=[0],
        z=[0],
        showlegend=False,
        mode='markers',
        marker=dict( size=3, color='grey', line={'width': 1, 'color': 'black'} )
        ))
    title = 'family: {}, mkey: {}, pkey: {}, ngj: {:.0f}'.format(
                                fam,
                                mkey,
                                data.id2meta[cellid]['pkey'],
                                data.gj_count[cellid])
    return {
    'data': traces,
    'layout': go.Layout(
        xaxis={'range': [-150, 150], 'zeroline': False, 'title':title},
        yaxis={'range': [-150, 150], 'zeroline': False},
        title=title,
        margin={'l': 40, 'b': 20, 't': 40, 'r': 10},
        height=height_px/2160*600
    )}  



# somatic Vm (of cells with DB) -----------------------------------------------
@app.callback(
    dash.dependencies.Output('graph-voltage', 'figure'),
    [dash.dependencies.Input('show-voltage', 'value'),
     dash.dependencies.Input('cellid', 'value')]
    )
def update_morph_graph(showVolt, cellid):
    
    if not showVolt: 
        return {}
    if not cellid:
        return {}
    v = data.voltage[cellid][::2]
    
    trace = [go.Scatter(
        x=data.time,
        y=v,
        text=cellid,
        mode='lines',
        opacity=0.5,
        showlegend=True,
        legendgroup=cellid,
        marker={ 'color': 'grey'},
        name=str(cellid)
        )]
    
    
    # depolarization block(s)?
    if cellid in data.cells_with_db:
        item = data.cells_with_db[cellid]
        x = []
        y = []
        for i,val in enumerate(item['start_indx']):
            e = data.time[item['end_indx'][i]]
            m = np.max(item['v'][val:item['end_indx'][i]])
            x += [data.time[val],e]
            y += [m,m]
            x.append(None)
            y.append(None)
        
        trace.append(go.Scatter(
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
    
                
    return {
            'data': trace,
            'layout': go.Layout(
                xaxis={'type': 'linear', 'title': 'time'},
                yaxis={'title': 'Vm'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 0},
                hovermode='closest'
                )
            }
    


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
      

    
if __name__ == '__main__':
    app.run_server(debug=True)
