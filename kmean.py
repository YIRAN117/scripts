#usage:
#python kmean.py <input Bacteria_Archaea file> <k value> <iteration time>
#example:python /kmean.py Bacteria_Archaea.tab 5 300
#output:
#WCSS,AIC,BIC,iteration number and group assigned to each species
import sys
import random
from math import *

def distEclud(vecA, vecB):
	return sqrt(sum([(vecA[i]-vecB[i])**2 for i in range(len(vecA))]))

def randcent(mat, k):
	ncol =len(mat[0])
	nrow=len(mat)
	ranges=[]
	centroid=[[0 for i in range(ncol)] for j in range(k)]
	for c in range(ncol):
		ranges.append([min([mat[i][c] for i in range(nrow)]),max([mat[i][c] for i in range(nrow)])])
	for n in range(k):
		for i in range(ncol):
			centroid[n][i]=ranges[i][0]+random.random()*(ranges[i][1]-ranges[i][0])
	return centroid

def kmeans(mat, k):
	ncol =len(mat[0])
	nrow=len(mat)
	clusterassign = [None for i in range(nrow)]
	d = [None for i in range(nrow)]
	centroids = randcent(mat, k)
	converge = False
	WCSS=0
	while not converge:
		converge = True
		centroid_new=[[0 for i in range(ncol)] for j in range(k)]
		cluster_count=[0]*k
		for r in range(nrow):
			min_dist = inf
			cluster_index = -1
			for n in range(k):
				dist = distEclud(centroids[n],mat[r])
				if dist < min_dist:
					min_dist = dist
					cluster_index = n
			if clusterassign[r] != cluster_index:
					converge = False
			cluster_count[cluster_index]+=1
			clusterassign[r] = cluster_index
			d[r]=min_dist
			for j in range(ncol):#update centroids
				centroid_new[clusterassign[r]][j]+=mat[r][j]
		for n in range(k):
			for i in range(ncol):
				if cluster_count[n]>0:
					centroids[n][i]=centroid_new[n][i]/cluster_count[n]
	for i in range(nrow):
		WCSS+=d[i]**2
	return centroids, clusterassign, WCSS

def kmean_best(mat,k,iter_max):
	ncol =len(mat[0])
	nrow=len(mat)
	WCSS_min = inf
	centroids_best = None
	clusterassign_best=None
	for i in range(iter_max):
	    centroids, clusterassign, WCSS = kmeans(mat,k)
	    if WCSS < WCSS_min:
	        WCSS_min = WCSS
	        centroids_best = centroids
	        clusterassign_best = clusterassign
	        iter_count=i
	bic=log(nrow)*k*ncol+WCSS_min
	aic=2*k*ncol+WCSS_min
	return  WCSS_min, aic, bic, centroids_best, clusterassign_best, iter_count

def normalize(mat):
	ncol =len(mat[0])
	nrow=len(mat)
	mat_n=centroid=[[0 for i in range(ncol)] for j in range(nrow)]
	s=[0]*ncol
	u=[0]*ncol
	sd=[0]*ncol
	ssm=[0]*ncol
	for i in range(ncol):
		for j in range(nrow):
			s[i]+=mat[j][i]
		u[i]=s[i]/nrow
	for i in range(ncol):
		for j in range(nrow):
			ssm[i]+=(mat[j][i]-u[i])**2
		sd[i]=sqrt(ssm[i]/nrow)
	for i in range(ncol):
		for j in range(nrow):
			mat_n[j][i]=(mat[j][i]-u[i])/sd[i]
	return mat_n

def main():
	mat=[]
	input_file=sys.argv[1]
	k=int(sys.argv[2])
	iter_max=int(sys.argv[3])
	with open(input_file,"r") as f1:
		for f in f1:
			if not f.startswith("Genome"):
				line=f.rstrip().split("\t")
				features=list(map(float, line[1:]))
				mat.append(features)
	mat_n=normalize(mat)
	ncol =len(mat_n[0])
	nrow=len(mat_n)
	WCSS_min, aic, bic, centroids_best, clusterassign_best,iter_count=kmean_best(mat_n,k,iter_max)
	print('{}\t{}\t{}\t{}\t{}'.format(WCSS_min, aic, bic,iter_count,clusterassign_best))
main()

