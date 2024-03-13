word MembrVolmodel(electrolyteProperties.lookup("MembraneVolumeModelling"));
if (MembrVolmodel == "yes")
{
    scalar MemVol(0);
    scalar N(0); 
    forAll(mesh.cellZones(), i)
    {
        const labelList& cellLabels = mesh.cellZones()[i];
        forAll(cellLabels, i)
        {
            label cellID = cellLabels[i]; 
            MemVol += mesh.V()[cellID];
        }
        scalar wsc_ = wsc.value();
        N = wsc_/MemVol;
        forAll(cellLabels, i)
        {
            label cellID = cellLabels[i];
            blockM.source()[cellID](2) += -N;
        }
    }
}