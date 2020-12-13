import pandas as pd
import numpy as np
import pickle
import plotly
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from scipy import signal

def visulize():
	with open('result_50familiy.data', 'rb') as filehandle:
		result_list = pickle.load(filehandle)

	df = pd.DataFrame(result_list, columns=["num_of_family", "num_of_structures", "num_of_sequences", "score_list"])
	df = df.drop(['score_list'], axis=1)
	# print(df)

	list = []
	for i in range(len(result_list)):
		a = pd.DataFrame(np.array(result_list[i][3]).reshape(1, len(result_list[i][3])))
		list.append(a)
	score = pd.concat([i.reset_index(drop=True) for i in list], axis=0)
	score.columns=["birch", "agglomerative", "k_means", "miniBatchKMeans", "dbscan", "meanShift", "spectralClustering", "gaussianMixture"]

	# print(score)

	df = pd.concat([df.reset_index(drop=True), score.reset_index(drop=True)], axis=1)
	print(df)
	bar_average_cluster(df)
	line_single_num_families(df)


def line_single_num_families(df):
	df = df.copy()
	fig = px.line(df, y="k_means", x="num_of_sequences", color="num_of_structures", title="k_means")
	# fig.write_image("line_single_num4_families.pdf")
	fig.show()
 
	df = df.copy()
	fig = px.line(df, y="agglomerative", x="num_of_sequences", color="num_of_structures", title="agglomerative")
	# fig.write_image("line_single_num4_families.pdf")
	fig.show()
 
	# df_4_family = df.loc[df['num_of_family'] ==4]
	# df_5_family = df.loc[df['num_of_family'] ==5]
	# df_6_family = df.loc[df['num_of_family'] ==6]

	# fig = px.line(df_4_family, y="k_means", x="num_of_sequences", color="num_of_structures", title="k_means #family=4")
	# fig.write_image("line_single_num4_families.pdf")
	# fig.show()

	# fig = px.line(df_5_family, y="k_means", x="num_of_sequences", color="num_of_structures", title="k_means #family=5")
	# fig.write_image("line_single_num5_families.pdf")
	# fig.show()

	# fig = px.line(df_6_family, y="k_means", x="num_of_sequences", color="num_of_structures", title="k_means #family=6")
	# fig.write_image("line_single_num6_families.pdf")
	# fig.show()


def bar_average_cluster(df):
	df = df.copy()
	list = []
	list.append(np.array(df['birch']).mean())
	list.append(np.array(df['agglomerative']).mean())
	list.append(np.array(df['k_means']).mean())
	list.append(np.array(df['miniBatchKMeans']).mean())
	list.append(np.array(df['dbscan']).mean())
	list.append(np.array(df['meanShift']).mean())
	list.append(np.array(df['spectralClustering']).mean())
	list.append(np.array(df['gaussianMixture']).mean())
	data = pd.DataFrame(list, columns=['score'])

	data['cluster']=["birch", "agglomerative", "k_means", "miniBatchKMeans", "dbscan", "meanShift", "spectralClustering", "gaussianMixture"]
	
	fig = px.bar(data, y='score', x='cluster', color='score', color_continuous_scale='Rainbow', title='Average for all clusters')
	fig.write_image("bar_average_cluster.pdf")
	fig.show()

	
