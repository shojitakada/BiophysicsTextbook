import numpy as np
import matplotlib.pyplot as plt

# parameter
b_transcript = 1.
b_translate = 0.3
a_m = 0.2
a_p = 0.04
maxiter = 10000

# initial condition
time = 0.
mRNA = 0
protein = 0
tstay = 0.
data = np.zeros((4, maxiter+1))
data[:,0] = [time, tstay, mRNA, protein]

# main loop
for iter in range(maxiter):
    rand1 = np.random.rand()
    rand2 = np.random.rand()
    option = []    # store expected changes in mRNA and protein in each option
    probability = []    # relative probability of each option


    op_transcript = [mRNA+1,protein] # transcription
    option.append(op_transcript)
    probability.append(b_transcript)

    op_translate = [mRNA,protein+1] # translate
    option.append(op_translate)
    probability.append(mRNA*b_translate)

    op_decay_m = [mRNA-1,protein] # decay of mRNA
    option.append(op_decay_m)
    probability.append(mRNA*a_m)

    op_decay_p = [mRNA,protein-1] # decay of protein
    option.append(op_decay_p)
    probability.append(protein*a_p)

    #print("probability", probability)

    tstay = -np.log(rand1)/np.sum(probability)    # propagation time
    time += tstay    # accumulated time

    prob_relative = probability/np.sum(probability)
    cum=0
    for k in range(len(probability)):
        cum+=prob_relative[k]
        if rand2 <= cum:
            winner = option[k]
            break
    winner = np.append([time,tstay],winner)
    mRNA   = winner[2]
    protein= winner[3]
    
    #print("winner",winner)  # for debug
    
    data[:,iter] = winner

    if time >= 10000:
        break


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data[0,0:iter], data[2,0:iter], color="b",  label="mRNA")
ax1.set_title("Number of mRNAs")
#fig1.savefig("mRNA.png")
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(data[0,0:iter], data[3,0:iter], color="r",  label="protein")
ax2.set_title("Number of proteins")
#fig2.savefig("protein.png")
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.scatter(data[2,np.int(iter/2):iter], data[3,np.int(iter/2):iter])
ax3.set_title("Correlataion of mRNA and protein")
#fig3.savefig("colleration.png")

