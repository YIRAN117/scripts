import sys
sys.setrecursionlimit(10000)
##compare if there is a match
def compare(a,b,match_score,mismatch_score):
	if a==b:
		value=match_score
	else:
		value=mismatch_score
	return value


##assign score based on the above,left and diagonal score 
##while store the position from where the new score inferred 
def assign(mat,i,j,match_score,mismatch_score,gap_score):
	y=mat[i-1][j]+gap_score
	z=mat[i][j-1]+gap_score
	w=mat[i-1][j-1]+compare(mat[i][0],mat[0][j],match_score,mismatch_score)
	value=y
	n=[i-1,j]
	if z>value:
		value=z
		n=[i,j-1]
	if w>value:
		value=w
		n=[i-1,j-1]
	return value,n

##trace back after matix filled
def path_back(mat,i,j,x_align,y_align,match_score,mismatch_score,gap_score):
	n=[i,j]
	if n[0]<2 or n[1]<2:
		return x_align,y_align
	else:
		n=assign(mat,i,j,match_score,mismatch_score,gap_score)[1]
		if n[0]==i-1 and n[1]==j-1:
			x_value=mat[0][j]
			y_value=mat[i][0]
		elif n[0]==i and n[1]==j-1:
			x_value=mat[0][j]
			y_value="-"
		else:
			x_value="-"
			y_value=mat[i][0]
		x_align.append(x_value)
		y_align.append(y_value)
		path_back(mat,n[0],n[1],x_align,y_align,match_score,mismatch_score,gap_score)

##print format
def insert_newlines(string, every=60):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

##read in files
x_file=sys.argv[1]
y_file=sys.argv[2]
match_score=float(sys.argv[3])
mismatch_score=float(sys.argv[4])
gap_score=float((sys.argv[5]))
with open(x_file,"r") as x:
	x_name=x.readline().rstrip().replace(">","")
	x_seq=x.read().replace("\n","").upper()
with open(y_file,"r") as y:
	y_name=y.readline().rstrip().replace(">","")
	y_seq=y.read().replace("\n","").upper()
x_len=len(x_seq)
y_len=len(y_seq)

##build the initial matrix
mat=[[None for i in range(x_len+2)] for j in range(y_len+2)]
mat[1][1]=0
for i in range(x_len):
	mat[0][i+2]=x_seq[i]
	mat[1][i+2]=mat[1][i+1]+gap_score
for j in range(y_len):
	mat[j+2][0]=y_seq[j]
	mat[j+2][1]=mat[j+1][1]+gap_score

##fill the matrix
for i in range(2,y_len+2):
	for j in range(2,x_len+2):
		mat[i][j]=assign(mat,i,j,match_score,mismatch_score,gap_score)[0]

##track back
x_align=[]
y_align=[]
path_back(mat,y_len+1,x_len+1,x_align,y_align,match_score,mismatch_score,gap_score)

##polish alignment ends
if len("".join(x_align).replace("-",""))<x_len:
	for i in reversed(range(0,x_len-len("".join(x_align).replace("-",""))-1)):
		x_align.append(x_seq[i])
		y_align.append('-')
if len("".join(y_align).replace("-",""))<y_len:
	for i in reversed(range(0,y_len-len("".join(y_align).replace("-",""))-1)):
		y_align.append(y_seq[i])
		x_align.append('-')
x_align=x_align[::-1]
y_align=y_align[::-1]
x_align_len=len(x_align)
y_align_len=len(y_align)

##create alignment quality line
middle_line=[None for i in range(x_align_len)]
for i in range(0,x_align_len):
	if x_align[i]==y_align[i]:
		middle_line[i]='*'
	else:
		middle_line[i]=' '

##print format
print("Sequence 1\t({}):\n{}\nlength:\t{}\n".format(x_name,insert_newlines(x_seq),x_len))
print("Sequence 2\t({}):\n{}\nlength:\t{}\n".format(y_name,insert_newlines(y_seq),y_len))
print("Scores:\nMatch: {}, Mismatch: {}, Gap: {}\n".format(match_score,mismatch_score,gap_score))
print("Alignment score: {}\nAlignment:\n".format(mat[-1][-1]))
a=1
for i in range(0, x_align_len, 60):
	space=' '*(len(str(a))+1)
	print('{}:{}\n{}{}\n{}:{}\n\n'.format(a,''.join(x_align[i:i+60]),space,''.join(middle_line[i:i+60]),a,''.join(y_align[i:i+60])))
	a+=60

