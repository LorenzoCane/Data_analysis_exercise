import numpy as np 
import matplotlib.pylab as plt
import glob

Use_Pixel_Mask = False
#Average = np.loadtxt('Threshold_prova.txt').reshape(-1,48,48)
Average = 1.05

#filename = 'inputpy10_bg1.txt'  #'eff_output_bg1prova_2ev.txt' # Path and filename of the txt file to be read
											# It requires a txt file with an integer number of 128 GTUs, 48x48 matrix
											# It skips the rows starting with 'Event'. That can be changed, A new event is expected to start after 128 GTUs


def compute_threshold__32gtu(D1, Average='32_gtu', Use_Pixel_Mask=False ):
	"""
	Set the threshold of each pixel 16 sigma above the average value in the packet. The average is computed on the first 32 GTUs of each packet.
	Returns the values of the threshold for each pixel in each packet. Threshold is made of K 48K48 matrix, where K is the number of D1 packets.
	"""
	if (type(Average)==str and Average=='32_gtu'):          
		print('using the first 32 GTUs to estimate the background')
		avg = np.mean(D1[:,0:32,:,:], axis=1)

	elif (type(Average)==int or type(Average)==float):          
		print('using ' + str(Average) + ' as background for all the pixels')
		avg = Average*np.ones((D1.shape[0],48,48))
		
	elif Average.shape[1:]==(48, 48):
		print('using  the values in the matrix Average as background for all the pixels')
		avg = Average

	else:
		print('Check the variable "Average"')
		print('It can be a number, a [nx48x48] matrix or the string "32_gtu"')

	threshold__32gtu = np.trunc(8*avg+(4*np.trunc(np.sqrt(128*avg)))) 
	is_below_28__32gtu = np.less(threshold__32gtu,28) 
	threshold__32gtu[is_below_28__32gtu] = 28

	if Use_Pixel_Mask == True:
		print('Using Pixel mask: All border pixels and the two bitshifted PMT')
		Pixel_mask = np.zeros((48,48),int)
		border = [0,7,8,15,16,23,24,31,32,39,40,47]
		PMT_11=np.array([[40,48],[32,40]]) 
		PMT_25=np.array([[8,16],[8,16]]) 
		for l in border:
			Pixel_mask[l,:]=1 
			Pixel_mask[:,l]=1 
		Pixel_mask[PMT_11[0,0]:PMT_11[0,1],PMT_11[1,0]:PMT_11[1,1]] = 1 
		Pixel_mask[PMT_25[0,0]:PMT_25[0,1],PMT_25[1,0]:PMT_25[1,1]] = 1 

		for packet in range (threshold__32gtu.shape[0]):
			for x in range (48):
				for y in range (48):
					if (Pixel_mask[x,y]==1):
						threshold__32gtu[packet,x,y] = 9999

	return threshold__32gtu


filenamelist = []
output_filenamelist = []

with open("inPY_bg1.txt","r") as tf:
	lines = tf.read().split('\n')
for line in lines:
	filenamelist.append(line)
	print(line)

for name in filenamelist:
	output_filenamelist.append('output_'+name) 

for file_number in range(len(filenamelist)-1):	
	f = open(filenamelist[file_number])     # open and put in a variable. The variable is an object 
	_data = f.read()
	f.close()       

	Events = _data.split('Event ')[1:]
	Data = np.zeros(len(Events), dtype=object) #visto che gli eventi hanno lunghezze diverse non posso metterli dentro una matrice normale. Quindi creo un vettore di oggetti, ogni oggetto nel vettore è un evento 

	n_triggers = np.zeros(len(Events)) # variable to store the number of pixels over threshold
	n_Event = np.zeros(len(Events)) # variable to store the event number
	for i in range(len(Events)): 
		n_Event[i] = int((Events[i][0:2])) 
		Data[i] = np.fromstring(Events[i].replace("\n", " "), dtype=int, sep=' ')[1:].reshape(-1,128,48,48)

		Threshold = compute_threshold__32gtu(Data[i], Average, Use_Pixel_Mask = False) # The value of the threshold of each pixel

		Boolean_matrix = np.empty((Data[i].shape[0],120,48,48), bool) # 1 if a pixel is over threshold, 0 otherwise 
		Sum8 = np.zeros((Data[i].shape[0],120,48,48))
		for gtu in range (120):
			Sum8[:,gtu,:,:] = Data[i][:,gtu,:,:]+Data[i][:,gtu+1,:,:]+Data[i][:,gtu+2,:,:]+Data[i][:,gtu+3,:,:]+Data[i][:,gtu+4,:,:]+Data[i][:,gtu+5,:,:]+Data[i][:,gtu+6,:,:]+Data[i][:,gtu+7,:,:].astype(np.uint16)
			Boolean_matrix[:,gtu,:,:] = np.greater(Sum8[:,gtu,:,:],Threshold)
		Matrix = np.sum(Boolean_matrix,axis=1) #How many times a pixels has been over threshold in each event
		n_triggers[i] = int(np.sum(Matrix))
	#per sapere quale pixel
	#np.where(Boolean_matrix==True)
	#PER SALVARE
	#np.savetxt("nome.txt",Boolean_matrix)

		GTU_in_packet = np.sum(Boolean_matrix, axis=(2,3)).astype(np.uint16) 
		Number_of_riggered_events = np.sum(np.sum(Boolean_matrix, axis=(1,2,3))!=0) 

		if (np.sum(Matrix)>0): #print only if there is a trigger
			print('In event ' + str(int(n_Event[i])) + ' in file ' + str(int(file_number) ) + ' there are ' + str( np.sum(Matrix)) + ' triggers \n')


	#Output_array = np.stack((n_Event,n_triggers),axis=1).astype(int) #Array saved into the txt output file. The first column is the number of event, the second the number of pixels above threshold
	#print('TEST ',n_triggers)
	#print('TEST ',len(n_triggers))
	print(file_number)
	arrayFileNum = np.full(len(n_Event),file_number)
	Output_array = np.stack((arrayFileNum,n_Event,n_triggers),axis=1).astype(int)
#	print('first ',Output_array)
#	np.insert(Output_array,0,file_number)
#	print('second ',Output_array)
#	Output_array.insert(0,file_number)
	print(file_number)
	np.savetxt(output_filenamelist[file_number], Output_array, fmt='%d')
	#print('TEST ',triggers)
	#np.savetxt(output_filenamelist[file_number],n_triggers, fmt='%d')



