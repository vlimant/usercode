def customise(process):
    from Workspace.MuonHLTTreeUtility.muonHLTTreeUtilityForData_cff import insertMHTU
    insertMHTU(process)
    return (process)
