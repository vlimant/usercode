def customise(process):
    from Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cff import insertMHTU
    insertMHTU(process)
    return (process)
