#!/bin/bash

reldir=$CMSSW_BASE/src
workdir=$reldir

#SOURCElist="source_adam.cff source_ttb.cff source_bjet.cff source_mbias.cff source_smu.cff source_4mu.flat.cff source_4mu.expo.cff"
#MAXEV=-1
#MODElist="fakeSTA.cff STA.cff correctedSTA.cff replacedSTA.cff reSTA.cff correctedreSTA.cff replacedreSTA.cff"

#SOURCElist="source_adam.cff source_ttb.cff source_bjet.cff source_mbias.cff"
#MAXEV=-1
#MODElist="fakeSTA.cff STA.cff correctedSTA.cff replacedSTA.cff"

#SOURCElist="source_4mu.flat.cff source_4mu.expo.cff source_adam.cff source_ttb.cff source_bjet.cff source_mbias.cff"
#MAXEV=50000
#MODElist="fakeSTA.cff reSTA.cff correctedreSTA.cff replacedreSTA.cff"
#MODElist="reSTA.cff correctedreSTA_F1.cff correctedreSTA_F2.cff"
#MODElist="reL2.cff correctedreL2_F1.cff correctedreL2_F2.cff"

SOURCElist="source_4mu.flat.cff source_muPlu10.180p1.cff"
MAXEV=-1
MODElist="correctedreL2_F1L2.cff reL2.cff correctedreOnFly_F1L2.cff correctedreOnFly_IdL2.cff"
TRACKERlist="localTracker_faster.cff localTracker_ondemand.cff localTracker_regular.cff"

for MODE in $MODElist ; do
for SOURCE in $SOURCElist ; do
for TRACKER in $TRACKERlist ; do

    dir=$workdir/$SOURCE.dir/$MODE.$TRACKER.dir
    mkdirhier $dir
    script=$dir/run.sh
    cfg=$dir/run.cfg

    ROOTOUTPUT=L3Muon.root

    cp $reldir/Configuration/L3Muon/data/L3MuonProducer.cfg.def $cfg
    replace -s MAXEV $MAXEV -- $cfg
    replace -s SOURCE $SOURCE -- $cfg
    replace -s ROOTOUTPUT $ROOTOUTPUT -- $cfg
    replace -s MODE $MODE -- $cfg
    replace -s TRACKER $TRACKER -- $cfg

echo $cfg

    cat > $script <<EOF

#!/bin/bash

#go to private code environement
cd ${reldir}
eval \`scramv1 runtime -sh\`

#large /pool scratch area
cd -

cmsRun ${dir}/run.cfg

#do timing
\$CMSSW_RELEASE_BASE/test/\$SCRAM_ARCH/analyzeTiming ${ROOTOUTPUT} L3Muon 1000 10000

echo ---------------------------------
echo printing the end of the log file
echo ---------------------------------
tail -200 detailedInfo.txt 
echo
echo last record seens
grep record detailedInfo.txt | tail -10
echo
echo

echo "in case you forgot something."
ls -trlh
echo "copying things to ${dir}/"
echo
#copy back things
cp HLT_*.root ${dir}/.
cp DQM_*.root ${dir}/.
cp *_after.root ${dir}/.
cp *_before.root ${dir}/.

/afs/cern.ch/user/v/vlimant/scripts/status.sh detailedInfo.txt

EOF
chmod 755 $script

done
done
done
