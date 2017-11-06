# This is a function which plots a comparison of the FDM and the ROM.
# It should plot side by side, titles, axis titles, and in the correct typeface.
def printing( model1=X,model1Title='High Fidelity Model',model2=U_ROM,model2Title='Data Driven Reduced Order Model',nLines=4 ):
	plt.subplot( 2,1,1 )
	plt.title( model1Title )
	plt.xtitle( 'Position' )
	plt.ytitle( 'Amplitude' )
		plt.plot( positionVector, model1[line*int(tstep/nLines)] )
	plt.subplot( 2,1,2 )
	plt.title( model2Title )
	plt.xtitle( 'Position' )
	plt.ytitle( 'Amplitude' )
		plt.plot( positionVector, model2[line*int(tstep/nLines)] )
	plt.show()
