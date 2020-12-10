#pragma rtGlobals=1		// Use modern global access method.
// Calcium analysis program by Thomas Voets (2011) - version 4.3 September 2017
// Use correct calibration constants!!
// WARNING: Routines in this file can kill wave files.
//	See the warnings below. Make sure you understand them before using these procedures.
//	NOTE: These functions works on waves in the current data folder only.

#pragma rtGlobals=1

#include <Strings as Lists>
#include <Execute Cmd On List>

// KillMatchingWaves(matchStr, flags)
//	matchStr is "*" for all waves or "wave*" for all waves whose name starts with "wave".
//	flags is:
//		bit 0: set to kill wave files (KillWaves/F)
//	WARNING: This can kill wave files. Make sure you know what you are doing.
//	NOTE: This function works on waves in the current data folder only.
Function KillMatchingWaves(matchStr, flags)
	String matchStr
	Variable flags
	
	String cmdTemplate
	
	if (flags %& 1)
		cmdTemplate = "KillWaves/Z/F %s"
	else
		cmdTemplate = "KillWaves/Z %s"
	endif
	ExecuteCmdOnQuotedList(cmdTemplate, WaveList(matchStr, ";", ""))
End

// RemoveWavesFromGraph(graphName, matchStr)
//	This is a building-block for RemoveWavesFromWindow.
//	NOTE: This function will not work reliably if the graph contains waves from other than the current data folder.
Function RemoveWavesFromGraph(graphName, matchStr)	
	String graphName						// name of a graph
	String matchStr						// "*" to remove all waves
	
	String wn
	String wl
	Variable i, offset, type
	
	DoWindow/F $graphName
	if (V_flag)
		wl = WaveList(matchStr, ";", "WIN:")			// list of waves in graph and in current data folder
		i = 0
		do
			wn = WaveName(graphName,i,1)			// next Y wave in graph
			if (strlen(wn) == 0)
				break										// all done
			endif
			offset = FindItemInList(wn, wl, ";", 0)
			if (offset >=0)
				RemoveFromGraph $wn
			else
				i += 1
			endif
		while (1)
	endif
End

// RemoveWavesFromTable(tableName, matchStr)
//	This is a building-block for RemoveWavesFromWindow.
//	NOTE: This function will not work reliably if the table contains waves from other than the current data folder.
Function RemoveWavesFromTable(tableName, matchStr)
	String tableName						// name of a table
	String matchStr						// "*" to remove all waves
	
	String wn
	String wl
	Variable i, offset
	
	DoWindow/F $tableName
	if (V_flag)
		wl = WaveList(matchStr, ";", "WIN:")			// list of waves in table
		i = 0
		do
			wn = WaveName(tableName,i,3)				// name of next column in table
			if (strlen(wn) == 0)
				break										// all done
			endif
			wn = wn[0, strsearch(wn, ".", 0)-1]		// get rid of suffix (.x, .y)
			if (CmpStr(wn[0], "'") == 0)				// remove single quotes if necessary
				wn = wn[1, strlen(wn)-2]
			endif
			offset = FindItemInList(wn, wl, ";", 0)
			if (offset >=0)
				RemoveFromTable $wn.id
			else
				i += 1
			endif
		while (1)
	endif
End

// RemoveWavesFromWindow(windowName, matchStr)
//	Removes waves whose names match the matchStr parameter from the named graph or table window.
//	NOTE: This function will not work reliably if the window contains waves from other than the current data folder.
Function RemoveWavesFromWindow(windowName, matchStr)
	String windowName
	String matchStr						// "*" to remove all waves
	
	DoWindow/F $windowName
	if (V_flag)
		if (WinType(windowName) == 1)			// this is a graph?
			RemoveWavesFromGraph(windowName, matchStr)
		endif
		if (WinType(windowName) == 2)			// this is a table?
			RemoveWavesFromTable(windowName, matchStr)
		endif
	endif
End

// RemoveWavesFromWindows(winMatchStr, matchStr)
//	Removes matching waves from matching windows.
//	winMatchStr is "*" for all waves or "Graph*" for all waves whose name starts with "Graph".
//	matchStr is "*" for all waves or "wave*" for all waves whose name starts with "wave".
//	NOTE: This function will not work reliably if the windows contain waves from other than the current data folder.
Function RemoveWavesFromWindows(winMatchStr, waveMatchStr)
	String winMatchStr					// "*" for all windows					
	String waveMatchStr					// "*" to remove all waves
	
	String cmdTemplate
	sprintf cmdTemplate, "RemoveWavesFromWindow(\"%%s\", \"%s\")", waveMatchStr
	ExecuteCmdOnList(cmdTemplate, WinList(winMatchStr,";","WIN:3"))
End

// KillAllWaves(flags)
//	Optionally, removes all waves from all graphs and tables.
//	Then kills all waves that are not in use.
//	flags is:
//		bit 0: set to kill wave files (KillWaves/F)
//		bit 1: set to remove all waves from graphs and tables
//	WARNING: This can kill wave files. Make sure you know what you are doing.
//	NOTE: This function works on waves in the current data folder only.
Function KillAllWaves(flags)
	Variable flags
	
	if (flags %& 2)
		RemoveWavesFromWindows("*", "*")
	endif
	KillMatchingWaves("*", flags)
End
	
	
	




Function MeanRow (wave2D, rn)
Wave wave2D
Variable rn			//row number

	MatrixOP/O/Free w=col(wave2D,rn)^t
	return mean(w)
End



function displayselecttraces(test): buttoncontrol
string test
WAVE todisplay=root:todisplay, responder=root:responder, nottodisplay=root:nottodisplay, problematic=root:problematic, problematictodisplay=root:problematictodisplay
WAVE GFPOK=root:GFPOK, badROI=root:badROI, basalcalciumOK=root:basalcalciumOK, relevant=root:relevant
string checkboxname, checkboxname2
variable numchecks= numpnts(events0)-1
variable ru=0
make/O/N=(numchecks) checkresults, checkresults2, checkresults3
variable count=1
todisplay=0
test =""
silent 1; pauseupdate
colortab2wave spectrumblack
do
checkboxname="check"+num2str(count)
checkboxname2="check2"+num2str(count)
controlinfo/W=selectionpanel $checkboxname
checkresults[count]=V_value
controlinfo/W=selectionpanel $checkboxname2
checkresults2[count]=V_value
controlinfo/W=selectionpanel removeugly
ru=V_value
count+=1
while (count<numchecks)
checkresults3=checkresults*checkresults2
duplicate/O responder, comparepos,compareneg, compareproblematic, compareboth, comparetotal
comparetotal=0
comparetotal=checkresults3[p]
comparepos[][]=(responder[p][q]>=checkresults[p])
todisplay[]=(meanrow(comparepos,p)==1)
compareneg[][]=(1-(checkresults2[p]*responder[p][q]==1))
nottodisplay[]=(meanrow(compareneg,p)==1)
compareboth=max(comparepos*compareneg,comparetotal)




compareproblematic[][]=max(problematic[p][q]*checkresults[p], problematic[p][q]*checkresults2[p])
problematictodisplay[]=1-(meanrow(compareproblematic,p)>0)


//todisplay *=nottodisplay

todisplay[]=(meanrow(compareboth,p)==1)
todisplay *=GFPOK
todisplay *=basalcalciumOK
todisplay *=badROI

relevant=GFPOK*basalcalciumOK*badROI
if(ru==1)
todisplay *=problematictodisplay
relevant *=problematictodisplay
endif




todisplay[0]=1


silent 1; pauseupdate
string name=""
variable teller=1 
variable /G ondisplay
Dowindow/K selection
Display/N=selection calcium0 vs nikontime 







do


name="calcium"+num2str(teller)
if(todisplay[teller]==1)
	AppendToGraph  $name vs nikontime
endif
	

teller+=1
while (exists("calcium"+num2str(teller)))
movewindow 0,40,600,350

NVAR  ishetratio=root:ishetratio


Label bottom "Time (minutes)"
if (ishetratio==0)
Label left "[Ca] (nM)"
SetAxis left 0,*

else
Label left "ratio"
SetAxis left 0,10
endif

DoWindow/T selection,"Selected calcium traces"
button/Z totable pos={300,20}, size={100,25}, proc=selecttracestottable, title="selection to table"


popupmenu whattotable pos={400,20}, value = "calcium traces;calcium amplitude;mean calcium;peak calcium;lowest calcium;maximal dCa/dt;mean dCa/dt;calcium SD;dCa/dt SD; End calcium; Number of Peaks" 





if (exists("events0"))
AppendToGraph/R events3 vs events0

ModifyGraph mode(events3)=1,rgb(events3)=(0,65280,65280)
AppendToGraph/R events3 vs events0
ModifyGraph mode(events3#1)=3; DelayUpdate
ModifyGraph lsize(events3)=0.2
ModifyGraph textMarker(events3#1)={events2,"default",1,90,5,0.00,-50.00}
ModifyGraph rgb(events3#1)=(0,0,0)
SetAxis right 0,1
ModifyGraph noLabel(right)=2
ModifyGraph tick(right)=3
endif


ondisplay = itemsinlist(wavelist("calcium*",";","WIN:selection"))
TextBox/C/N=text0/A=MC "Number of cells: "+num2str(ondisplay-1)+" out of "+num2str(sum(relevant)-1) +"  ("+num2str(round(100*(ondisplay-1)/(sum(relevant)-1)))+" %)"


modifygraph lSize=0.5

if(ondisplay>1)
teller=0
	do
	teller+=1
	while(todisplay[teller]==0)
	scrolldata("",teller,"calcium"+num2str(teller))
else 
dowindow/K scrolldisplay
endif

end

macro restore(test): buttoncontrol
string test

badROI=1
displayselecttraces("")
end



function findresponders (ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum	// value of variable as number
	String varStr		// value of variable as string
	String varName	// name of variable
	WAVE basalcalciumOK=root:basalcalciumOK, GFPOK=root:GFPOK, responder=root:responder, problematic=root:problematic, basalcalcium=root:basalcalcium, GFP=root:GFP, GFPcursorX=root:GFPcursorX
	WAVE maxdif=root:maxdif, meandif=root:meandif, sddif=root:sddif, amplitude=root:amplitude, maxcal=root:maxcal, endcal=root:endcal, sdcal=root:sdcal, prepeak=root:prepeak, fluodifproduct=root:fluodifproduct, numosc=root:numosc
	NVAR supercontrolcheckboxvalue, stringency, calciumthreshold, basalcallow, basalcalhigh, numwaves
	
variable wavecount, eventcount, supercontrol

basalcalciumOK=1
GFPOK=1
variable low2, high2
//low2 = min(low, high)
//high2= max(low,high)
low2=-65000
high2=65000

silent 1; pauseupdate
responder=0
problematic=0

supercontrol =1-supercontrolcheckboxvalue

wavecount=1
do
eventcount=1

	do
	
	if((maxdif[eventcount][wavecount]>(meandif[eventcount-1][wavecount]+stringency*sddif[0][wavecount]))&(amplitude[eventcount][wavecount]>calciumthreshold))
		responder[eventcount][wavecount]+=.4
	endif	
	if((maxcal[eventcount][wavecount]>(endcal[eventcount-1][wavecount]+stringency*sdcal[0][wavecount])))
		responder[eventcount][wavecount]+=.6
	endif
	if(prepeak[eventcount][wavecount]<0.0)
		responder[eventcount][wavecount] *= supercontrol
	endif
	
	if(fluodifproduct[eventcount][wavecount]<-.2)
	problematic[eventcount][wavecount] =1
	endif
	

	
	eventcount +=1
	while(eventcount<=numpnts(events0)-1)
	
if((basalcalcium[wavecount]<basalcallow) | (basalcalcium[wavecount]>basalcalhigh))
basalcalciumOK[wavecount]=0
endif

//print GFP[wavecount], GFPcursorX[0], GFPcursorX[1]
if((GFP[wavecount]<GFPcursorX[0]) | (GFP[wavecount]>GFPcursorX[1]))
GFPOK[wavecount]=0
endif

wavecount +=1
while (wavecount<=numwaves)	


displayselecttraces("")


end



macro  GFPexecute(test) : buttoncontrol
	String test
	

findresponders("",1,"","")

end

macro selecttracestottable(test): buttoncontrol
string test

pauseupdate
Make/O/N=(10,10)/T W_wavelist 
Getwindow selection, wavelist

string relevantwave

string ROIstring
variable wavenumber



Dowindow/K selectiontable


edit/N=selectiontable 

variable count=1

controlinfo whattotable
variable switchvalue = V_value
if (switchvalue==1)
appendtotable/W=selectiontable nikontime
endif
if (switchvalue>1)
appendtotable/W=selectiontable events0, events1, events2
endif



do
relevantwave=W_wavelist[count][0]
if(stringmatch(relevantwave,"event*")==0)
wavenumber=str2num((W_wavelist[count][0])[7,12])
	if(switchvalue==1)	
			ROIstring= "ROI_calcium"+num2str(ROI(wavenumber))
			duplicate/O  $W_wavelist [count][0], $ROIstring
			appendtotable/W=selectiontable $ROIstring
	endif					// exit from switch
	if(switchvalue==2)	
			ROIstring= "ROI_amplitude"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] amplitude, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
	endif
	if(switchvalue==3)		
			ROIstring= "ROI_meancal"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] meancal, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif
	if(switchvalue==4)			
			ROIstring= "ROI_maxcal"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] maxcal, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif			
	if(switchvalue==5)		
			ROIstring= "ROI_mincal"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] mincal, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif		
	if(switchvalue==6)		
			ROIstring= "ROI_maxdCadt"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] maxdif, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif		
	if(switchvalue==7)
			ROIstring= "ROI_meandCadt"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] meandif, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif		
	if(switchvalue==8)	
			ROIstring= "ROI_calciumSD"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] SDcal, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
			endif	
	if(switchvalue==9)		
			ROIstring= "ROI_dCadtSD"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] SDdif, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
		endif	
	
	if(switchvalue==10)		
			ROIstring= "ROI_endcal"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] endcal, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
		endif	
		
	if(switchvalue==11)		
			ROIstring= "ROI_peaks"+num2str(ROI(wavenumber))
			duplicate/O/R=[][wavenumber] numosc, $ROIstring
			redimension/N=(numpnts(events0)) $ROIstring
			appendtotable/W=selectiontable $ROIstring
		endif		

endif
count+=1
while(count<numpnts(W_wavelist))
movewindow/W=selectiontable 850,0,1420,355

resumeupdate
end



macro loadolympus340380(bgwave)
variable bgwave=1

XLLoadWave/R=(A3,ZZ10000)/Q/D/A=OF340 
killwaves/Z timeolympus
rename OF3400, timeolympus
timeolympus /=60*24
timeolympus[0]=0
concatenate/NP wavelist("OF340*", ";", ""), Intensity_340
duplicate/O Intensity_340, background_340

XLLoadWave/R=(A3,ZZ10000)/Q/D/A=OF380 
killwaves/Z OF3800
concatenate/NP wavelist("OF380*", ";", ""), Intensity_380
duplicate/O Intensity_380, background_380, ND_T, ROI_ID

string background340="OF340"+num2str(bgwave), background380="OF380"+num2str(bgwave)

background_380[] =$background380[mod(p, numpnts(timeolympus))]
background_340[] =$background340[mod(p, numpnts(timeolympus))]
ND_T[]=mod(p, numpnts(timeolympus))+1
ROI_ID[]=floor(p/numpnts(timeolympus))+1

background_380 *=10
Intensity_340 *=10
background_340 *=10
Intensity_380 *=10

end


macro calciumanalysis(Rmin, Rmax, alfa, filetype,events, bg340,bg380, ratioiocalcium)
variable Rmin = 0.323, Rmax = 10.59, alfa = 0.71 // 10x objective 
//variable Rmin = 0.71969, Rmax = 9.97446, alfa=0.8 // 20x objective
//variable Rmin = 0.61077, Rmax = 6.99908, alfa = ?? // 40x objective 
string filetype
prompt filetype, "Datafile", popup, "Olympusexcel;newNikon;Olympustxt"
string events
prompt events, "Eventfile", popup, "newNikon;none"
variable bg340=0
prompt bg340, "Enter value for F340 bg, or zero for ROI: "
variable bg380=0
prompt bg380, "Enter value for F380 bg, or zero for ROI: "
string ratioiocalcium
prompt ratioiocalcium, "What?", popup, "calcium;ratio"


print ratioiocalcium
variable Keff = 260*(Rmax+alfa)/(Rmin+alfa)

killallwaves(2)

silent 1; pauseupdate


variable/G backgroundwave=0,numwaves, backgroundroi=0
variable count=1, finalcount, numpoints,   backgroundgfpsd=0
string columnname, calciumname, calciumname2, difname, difsmoothname, difrationame, productname, f340wave, f380wave,f340bgwave,f380bgwave, fGFPwave, fGFPbgwave, timewave


if (cmpstr(filetype,"Olympusexcel")==0)

loadolympus340380(1)
else
XLLoadWave/Q/W=1/O
		

endif

//if (cmpstr(filetype,"newNikon")==0)

if ((1-1)==0)
f340wave = wavelist("Intensity*340*","","")
f380wave = wavelist("Intensity*380*","","")
fGFPwave= wavelist("Intensity*GFP*","","")
fGFPbgwave= wavelist("background*GFP*","","")
f340bgwave = wavelist("background*340*","","")
	if (strlen(f340bgwave)<2)
	abort "Background trace could not be found"
	endif

f380bgwave = wavelist("background*380*","","")
timewave  = wavelist("Time*","","")
duplicate/O $f340wave correctedratio, correctedf340, correctedf380, correctediso


	if (bg340==0)
	
		if ( (wavemin($f340wave)<wavemin($f340bgwave))||(wavemin($f380wave)<wavemin($f380bgwave)))
		beep
		doalert/T="Warning!" 1,"Fluorescence in background ROI is higher than in some other ROIs. Take lowest value as background?"
			if (V_flag==2)
			abort
			endif
		endif
	
	
	bg340=min(wavemin($f340bgwave), wavemin($f340wave))
	bg380=min(wavemin($f380bgwave), wavemin($f380wave))
	//print bg340, bg380
	//bg340=110
	//bg380=110
		
	endif
print "Used background (340,380)", bg340,bg380

correctedf340 = $f340wave -bg340
smooth/M=0 10,correctedf340
correctedf380 =$f380wave - bg380
smooth/M=0 10,correctedf380
correctedratio= correctedf340/correctedf380
correctediso = correctedf340+alfa*correctedf380
//correctedf340 /=correctediso
//correctedf380 /=correctediso

wavestats/Q ND_T
numpoints = V_max
numwaves = V_npnts/V_max

Make/D/O/N=(numwaves+1) GFP, ROI, badROI, noisyROI
GFP=0; badROI=1; noisyROI=0
count=1
	
	
	do
	ROI[count]=ROI_ID(numpoints*count-5)
	if (strlen(fGFPwave)>2)
	GFP(count)=$fGFPwave(numpoints*count-5)-$fGFPbgwave(numpoints*count-5)
	endif
	count+=1
	
	while (count<numwaves+2)
	



count=1
do
columnname = "trace"+num2str(count)
duplicate/R=[0+(count-1)* numpoints,count* numpoints-1]/O correctedratio $columnname


wavestats/Q $columnname
	if (V_numNaNs>1)
	backgroundwave=count
	backgroundroi=1
	print "backgroundwave = number: " + num2str(ROI(count)) + "calcium"+num2str(count)
		if (exists("Background_GFP_")==1)
		backgroundgfpsd=StDEV_GFP_(numpoints*count-5)
		endif
	endif
columnname = "f340trace"+num2str(count)
difname = columnname + "DIF"
productname = "product" + num2str(count)
duplicate/R=[0+(count-1)* numpoints,count* numpoints-1]/O correctedf340 $columnname
smooth/M=0 5,$columnname
Differentiate $columnname/D=$difname
duplicate $difname $productname

columnname = "f380trace"+num2str(count)
difname = columnname + "DIF"
duplicate/R=[0+(count-1)* numpoints,count* numpoints-1]/O correctedf380 $columnname
smooth/M=0 5,$columnname
Differentiate $columnname/D=$difname


$productname *=-$difname


		
columnname = "isotrace"+num2str(count)
duplicate/R=[0+(count-1)* numpoints,count* numpoints-1]/O correctediso $columnname

$productname /=mean($columnname)
integrate $productname

count +=1
while (count<numwaves+1)


duplicate/R=[0, numpoints-1]/O $timewave nikontime 

nikontime *=60*24

endif

if (cmpstr(events,"newnikon")==0)

XLLoadWave/Q/C=3/O/N=events
Make/N=(numpnts(events0))/D/O events3
count=0
do
	if(stringmatch(events1[count], "user*")==0)
	deletepoints count,1, events0, events1, events2, events3
	count -=1
	endif
	if((count>0)&(stringmatch(events1[count],events1[count-1])))
		deletepoints count,1, events0, events1, events2, events3
	count -=1
	endif
	
	count +=1
while (count<numpnts(events0))
events3 = 1
events0 /=60
	if(events0[0]>nikontime[1])
	InsertPoints 0,1, events0,events1,events2,events3
	events0[0]=nikontime[0]
	events1[0]="Userstart"
	events2[0]="start experiment"
	events3[0]=1
	else
	events0[0]=nikontime[0]
	endif
insertpoints (numpnts(events0)),1, events0,events1,events2,events3
wavestats/Q nikontime
events0[numpnts(events0)-1]=V_max-0.02
events1[numpnts(events0)-1]="Userend"
events2[numpnts(events0)-1]="end experiment"
events3[numpnts(events0)-1]=1

else //make events if no file is loaded
Make/N=2/T/O  events1, events2
Make/N=2/D/O events0, events3
events0[0]=nikontime[0]
events1[0]="Userstart"
events2[0]="start experiment"
events3[0]=1
wavestats/Q nikontime
events0[1]=V_max-0.0
events1[1]="Userend"
events2[1]="end experiment"
events3[1]=1


endif




count=1

do

columnname = "trace"+num2str(count)
calciumname="calcium"+num2str(count)
calciumname2="calciumb"+num2str(count)
difname="dif"+num2str(count)
difsmoothname ="difsmooth"+num2str(count)
difrationame="difratio"+num2str(count)




duplicate/O $columnname $calciumname, $calciumname2, tellercorrected, noemercorrected

tellercorrected=max(0,($columnname-Rmin))
tellercorrected=min(Rmax-Rmin,$columnname-Rmin)

noemercorrected=max((Rmax-$columnname),(Rmax-Rmin)*0.15)

$calciumname2 = Keff*($columnname-Rmin)/(Rmax-$columnname)
$calciumname = Keff*tellercorrected/noemercorrected
if(stringmatch(ratioiocalcium,"ratio"))
$calciumname2 = $columnname
$calciumname= $columnname
endif

$calciumname[0]=$calciumname[1]
duplicate/O $calciumname $difname

differentiate/METH=2 $calciumname/X=nikontime/D=$difname
duplicate/O $difname, $difsmoothname
smooth 10, $difsmoothname

$difname/=60
$difsmoothname/=60

duplicate/O $columnname $difrationame

differentiate/METH=2 $difrationame/X=nikontime

$difrationame/=60

count+=1
while (exists("trace"+num2str(count))==1)
duplicate/O calcium1, calcium0
calcium0=0

finalcount = numwaves

make/O/N=(finalcount) peaks, peakpos
count=1
do
difname="dif"+num2str(count)
wavestats/Q $difname
peaks[count]=V_max
peakpos[count]=V_maxloc

count+=1
while (count<finalcount)


Dowindow/K diffigure
Dowindow/K calciumfigure
Dowindow/K GFPhistogram
Dowindow/K selectionpanel






Make/N=200/O GFP_Hist;DelayUpdate
Histogram/B=1 GFP,GFP_Hist
dowindow /K GFPhistogram
Display /N=GFPhistogram GFP_Hist as "GFP histogram"
ModifyGraph mode=1,hbFill=2,rgb=(0,65280,0)
movewindow/W=GFPhistogram 610,430,1000,730
wavestats/Q GFP


//showinfo

variable/G low, high, GFPcount


Make/O/N=2 GFPcursorY, GFPcursorX
GFPcursorx[0]=V_min-(V_max-V_min)/10-1
GFPcursorx[1]=V_max+(V_max-V_min)/10+1
GFPcursorY=wavemax (GFP_hist)

appendtograph/W=GFPhistogram GFPcursorY vs GFPcursorX
ModifyGraph mode=1,lsize(GFPcursorY)=3

button/Z applyGFP win=GFPhistogram,pos={200,100},fcolor=(0,50000,0), size={160,25}, proc=GFPexecute, title="Apply GFP range"
TextBox/C/N=text0/F=0/A=MT/F=2 "Use Left and Right arrows to change maximum\rUse < and > to change minimum"

setwindow GFPhistogram, hook(myhook2)=mywindowhook2


//setdrawlayer/W=GFPhistrogram ProgBack
//setdrawenv xcoord=bottom, ycoord=left, linethick=2
//drawline low,0,low,20

//setdrawlayer UserBack
//setdrawenv xcoord=bottom, ycoord=left, linethick=2

//drawline high,0,high,20



newpanel/N=selectionpanel/W=(820,0,1120,240+numpnts(events1)*15) as "Select response profile"




button/Z displayit pos={40,10}, size={100,25}, proc=displayselecttraces, title="display selection" 



//slider testslider1 pos={50,0},size={440,50}, ticks=0, limits={V_min-100,V_max+100,10}, vert =0, variable =low, proc=setbounds1
//slider testslider2 pos={50,40},size={440,50}, ticks=0, limits={V_min-100,V_max+100,10},  vert =0, variable =high, proc =setbounds2




variable/G stringency=3, calciumthreshold=(0.05+49.95*stringmatch(ratioiocalcium,"calcium")), basalcallow=0, basalcalhigh=10000, supercontrolcheckboxvalue = stringmatch(ratioiocalcium,"calcium"), ishetratio= (1-stringmatch(ratioiocalcium,"calcium"))
count=1
string checkboxname, checkboxname2

string titlestring
do
checkboxname ="check"+num2str(count)
checkboxname2="check2"+num2str(count)
CheckBox $checkboxname,pos={100,150+15*count},size={78,15},title=events2[count],value= 0,mode=0
CheckBox $checkboxname2,pos={200,150+15*count},size={78,15}, title = "", value= 0,mode=0
count+=1
while (count<numpnts(events1)-1))
checkbox removeugly, pos={100,170+15*count},size={78,15},title="Disregard problematic?", value=stringmatch(ratioiocalcium,"calcium"), mode=0
button/Z bringback pos={75,195+15*count},fcolor=(50000,0,0), size={160,25}, proc=restore, title="bring back bad traces"

setvariable stringencyinput limits={1,75,.5}, pos={60,40},size={150,50}, title="Stringency",  proc = findresponders
setvariable stringencyinput value=stringency

titlestring = "Minimal increase"
if (stringmatch(ratioiocalcium,"calcium"))
titlestring +=" (nM)"
else
titlestring +=" (R)"
endif

setvariable calciumthresholdinput limits={0,10000,5}, pos={60,80},size={180,50}, title=titlestring,  proc = findresponders
setvariable calciumthresholdinput value=calciumthreshold

titlestring = "Basal"
if (stringmatch(ratioiocalcium,"calcium"))
titlestring +=" calcium"
else
titlestring +=" ratio"
endif


setvariable calciummininput limits={0,10000,5}, pos={60,100},size={140,50}, title=titlestring,  proc = findresponders
setvariable calciummininput value=basalcallow

setvariable calciumhighinput limits={0,10000,5}, pos={210,100},size={80,50}, title="and <",  proc = findresponders
setvariable calciumhighinput value=basalcalhigh

setvariable supercontrolcheckbox limits={0,1,1}, pos={60,60},size={100,50}, title= "supercontrol",   proc = findresponders
setvariable  supercontrolcheckbox value=supercontrolcheckboxvalue



Make/N=((numpnts(events0)-1), numwaves+1)/D/O meancal=0, maxcal=0, mincal=0, sdcal=0,begincal=0, endcal=0, amplitude=0, maxdif=0, meandif=0, numosc=0, sddif=0, responder=0, problematic=0, fluodifproduct=0, prepeak=0
Make/O/N=(numwaves+1) basalcalcium, qq, basalcalciumOK=1, GFPOK=1, relevant=1, todisplay=0, nottodisplay=0, problematictodisplay=0,anyhowtodisplay=0
variable wavecount=1, eventcount, rangemin, rangeplus


silent 1; pauseupdate
do
eventcount=0
calciumname = "calcium"+num2str(wavecount)
difname ="dif"+num2str(wavecount)
difsmoothname="difsmooth"+num2str(wavecount)

productname ="product"+num2str(wavecount)

	do
	findlevel/P/Q nikontime, events0[eventcount]
	rangemin=round(V_levelX)
	findlevel/P/Q nikontime, events0[eventcount+1]
	rangeplus = max (round(V_levelX)-1, rangemin + 7)
//print rangemin, rangeplus

	
	wavestats/Q/R=[rangemin, rangeplus]/Q $calciumname
	meancal[eventcount][wavecount]=V_Avg
	maxcal[eventcount][wavecount]=V_max
	
	
	mincal[eventcount][wavecount]=V_min
	sdcal [eventcount][wavecount]=V_SDev
		if(V_SDev<sdcal[0][wavecount])
		sdcal[0][wavecount]=V_SDev
		endif
	begincal[eventcount][wavecount]= mean ($calciumname,pnt2x($calciumname,rangemin), pnt2x($calciumname,rangemin+4))
	endcal[eventcount][wavecount]= mean ($calciumname, pnt2x($calciumname,rangeplus-4), pnt2x($calciumname,rangeplus))
	if(eventcount>0)
		amplitude[eventcount][wavecount]=	maxcal[eventcount][wavecount]-endcal[eventcount-1][wavecount]
	endif	
	wavestats/Q/R=[rangemin, rangeplus]/Q $difname
	maxdif[eventcount][wavecount]=V_max
	prepeak[eventcount][wavecount]=max($productname[V_maxrowloc]-$productname[V_maxrowloc-5], $productname[V_maxrowloc]-$productname[rangemin])
	sddif [eventcount][wavecount]=V_SDev
		if(V_SDev<sddif[0][wavecount])
		sddif[0][wavecount]=V_SDev
		endif
	meandif[eventcount][wavecount]=V_Avg
	
	
//Hier aanpassen als ge te veel of te weinig peaks vindt	
	//findlevels/B=5/EDGE=1/R=[rangemin-2,rangeplus-1]/Q $difsmoothname,max (0.3,2*sddif[0][wavecount])
	//numosc[eventcount][wavecount]=V_Levelsfound
	
	wavestats/Q/R=[rangemin, rangeplus]/Q $productname
	fluodifproduct[eventcount][wavecount]=($productname[rangeplus]-$productname[rangemin])

	
	
	eventcount +=1
	while(eventcount<=numpnts(events0)-2)
//meandif[0][wavecount]=mean($difname)
basalcalcium[wavecount]=meancal[0][wavecount]
qq[wavecount]=statsmedian($productname)

		if(qq[wavecount]<-1)
		noisyROI[wavecount]=1
		endif

wavecount +=1
while (wavecount<=numwaves)	


wavestats/Q/R=[1,] basalcalcium
basalcallow = floor(V_min)
basalcalhigh=ceil(V_max)
findresponders("",3,"","")

wavecount=0


end


function scrolldata(ctrlname,popNum, popStr): popupmenucontrol
string ctrlname
variable popNum
string popStr

Variable v1
sscanf popStr, "calcium %f", v1
string f340scroll="f340trace"+num2str(v1)
string f380scroll="f380trace"+num2str(v1)
string calciumscroll ="calcium"+num2str(v1)
string calciumscrollb ="calciumb"+num2str(v1)
dowindow/K scrolldisplay
variable/G scrollwavenumber=1
display/N=scrolldisplay $calciumscroll vs nikontime as "Scroll individual ROIs"

SetAxis left 0,*
appendtograph/R $f340scroll, $f380scroll vs nikontime

ModifyGraph/W=scrolldisplay tick(right)=0,noLabel(right)=0

movewindow 0,370,600,570

ModifyGraph lsize($calciumscroll)=2,rgb($calciumscroll)=(0,0,65280)
ModifyGraph lsize($f340scroll)=2,rgb($f340scroll)=(0,65280,33024)
ModifyGraph lsize($f380scroll)=2

Label bottom "Time (minutes)"
Label right "Fluorescence (a.u.)"

NVAR  ishetratio=root:ishetratio

if (ishetratio==0)
Label left "[Ca] (nM)"
textBox/C/N=text0/F=0/A=RT/y=-30 "\\Z16\\K(0,0,65280)Calcium\r\\K(0,65280,0)F340\r\\K(65280,0,0)F380"

else
Label left "ratio"
textBox/C/N=text0/F=0/A=RT/y=-30 "\\Z16\\K(0,0,65280)Ratio\r\\K(0,65280,0)F340\r\\K(65280,0,0)F380"

endif


if (exists("events0"))
AppendToGraph/R=events events3 vs events0

ModifyGraph mode(events3)=1,rgb(events3)=(0,65280,65280)
AppendToGraph/R=events events3 vs events0
ModifyGraph mode(events3#1)=3; DelayUpdate
ModifyGraph lsize(events3)=0.2
ModifyGraph textMarker(events3#1)={events2,"default",1,90,5,0.00,-50.00}
ModifyGraph rgb(events3#1)=(0,0,0)
SetAxis right 0,1
ModifyGraph noLabel(right)=2
ModifyGraph tick(right)=3
ModifyGraph freePos(events)=0
SetAxis events 0,1
ModifyGraph tick(events)=3
ModifyGraph axThick(events)=0
ModifyGraph noLabel(events)=2
endif
SetAxis right 0,*

ModifyGraph noLabel(right)=0,noLabel(events)=2
ModifyGraph/W=scrolldisplay margin(top)=60


popupmenu scrollwindow mode=1,pos={100,30}, title="Wave to display", popvalue=popStr, value = wavelist("calcium*",";","WIN:selection")[0,strlen(wavelist("calcium*",";","WIN:selection"))-10], proc=scrolldata
//button scrollup, pos ={270,20}, size={100,20}, title="Scroll Up", proc=scrollup
//button scrolldown, pos ={270,40}, size={100,20}, title="Scroll Down", proc=scrolldown
button/Z dbt pos={400,30}, size={100,25}, proc=deletebadtrace, title="delete bad trace"


setwindow scrolldisplay, hook(myhook)=mywindowhook

end

function scrollup(test): buttoncontrol
string test
Wave todisplay=root:todisplay

controlinfo scrollwindow

silent 1; pauseupdate

Variable v1,teller=0
v1=str2num(S_value[7,12])



teller=v1

do
	
teller+=1
if (teller==(numpnts(todisplay)))
teller=1
endif
while (todisplay[teller]==0)

scrolldata("",teller,"calcium"+num2str(teller))

end



function scrolldown(test): buttoncontrol
string test
Wave todisplay=root:todisplay

controlinfo scrollwindow
silent 1; pauseupdate
Variable v1,teller=0
v1=str2num(S_value[7,12])


teller=v1

do
	
teller-=1
if (teller==0)
teller=numpnts(todisplay)-1
endif
while (todisplay[teller]==0)

scrolldata("",teller,"calcium"+num2str(teller))

end




function deletebadtrace(test): buttoncontrol
string test
WAVE badROI = root:badROI

controlinfo scrollwindow

Variable v1,teller=0
v1=str2num(S_value[7,12])



badROI[v1]=0
displayselecttraces("")

end




function MyWindowHook(s)
	STRUCT WMWinHookStruct &s
	
	
	Variable hookResult = 0	// 0 if we do not handle event, 1 if we handle it.
	string/G etc="etc"
	switch(s.eventCode)
		case 11:					// Keyboard event
			switch (s.keycode)
				
				case 31:
					scrollup(etc)
					//print "up"
					hookResult = 1
					break
				case 30:
					scrolldown(etc)
					//print "down"
					hookResult = 1
					break			
			endswitch
			break
	endswitch

	return hookResult	// If non-zero, we handled event and Igor will ignore it.
End





function MyWindowHook2(s)
	STRUCT WMWinHookStruct &s
	WAVE GFPcursorX = root:GFPcursorX
	WAVE GFP = root:GFP
		
	Variable hookResult2 = 0	// 0 if we do not handle event, 1 if we handle it.
	string/G etc2="etc2"
	switch(s.eventCode)
		case 11:					// Keyboard event
		//print s.keycode
			switch (s.keycode)
				
				case 28:
					//scrollup(etc2)
					//print "rightleft"
					if(GFPcursorX[1]-wavemax(GFP)/50>GFPcursorX[0])
					GFPcursorX[1]-=wavemax(GFP)/50
					else
					GFPcursorX[1]=GFPcursorX[0]
					endif
					hookResult2 = 1
					break
				case 29:
					//scrolldown(etc2)
					//print "rightright"
					GFPcursorX[1]+=wavemax(GFP)/50
					hookResult2 = 1
					break			
				case 44:
					//scrollup(etc2)
					//print "leftleft"
					
					GFPcursorX[0]-=wavemax(GFP)/50
					
					hookResult2 = 1
					break
				case 46:
					//scrolldown(etc2)
					//print "rightright"
					if(GFPcursorX[0]+wavemax(GFP)/50<GFPcursorX[1])
					GFPcursorX[0]+=wavemax(GFP)/50
					else
					GFPcursorX[0]=GFPcursorX[1]
					endif
					
					
					hookResult2 = 1
					break			
		
					
					
					
			endswitch
			break
	endswitch

	return hookResult2	// If non-zero, we handled event and Igor will ignore it.
End