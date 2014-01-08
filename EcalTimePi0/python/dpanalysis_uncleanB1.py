import FWCore.ParameterSet.Config as cms


filelist = cms.untracked.vstring()

filelist.extend([
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/002F27BF-F371-E211-B732-003048FFCC2C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/00580670-9F71-E211-8180-003048678BB2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/007433CF-0072-E211-B710-00261894386F.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/007951A5-0272-E211-A13F-002618943834.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0079D246-F471-E211-ADD8-0026189438D6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/007A563F-F471-E211-9F26-00261894388D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/007D8B63-7071-E211-9D44-0025905938D4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/008BFBCD-FA71-E211-B473-003048678FD6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/00D0953F-BE71-E211-BB52-0025905938A4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/00D381DF-AC71-E211-A6A5-0026189438FF.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/00E61F3A-BA71-E211-BF7C-003048679266.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/00FBA5C4-FD71-E211-B1D6-002354EF3BE0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0209F0BD-F571-E211-89F4-0026189438D9.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/020CD57B-0472-E211-B24B-00261894389E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02201CFE-0872-E211-9018-0030486790A0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02206548-C171-E211-BCBC-002618943829.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0224A23B-0572-E211-9939-003048678B12.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0230954F-A271-E211-B7EB-003048678B0E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02548232-F471-E211-95ED-00261894394F.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/026463D4-7571-E211-88F2-0030486790B8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02648872-BC71-E211-B0FA-00261894387B.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0277DBD6-0572-E211-A70F-003048678BAC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/027A7E07-F371-E211-B8F3-002618943933.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02E0ACDA-FE71-E211-9F96-002618943849.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02F3FB82-E071-E211-9AE6-00259059649C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/02FC3C8E-DB71-E211-AF2A-003048FFD752.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/04049694-E071-E211-9162-0025905964CC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/040BC2CD-B571-E211-97F6-003048FFCB6A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/04135898-0172-E211-9C96-002618943882.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/042832B9-7971-E211-84AE-003048FFCB8C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0431E73C-C971-E211-A233-003048FFD770.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/04362B96-C771-E211-A385-003048FFD71A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/043B59E7-F271-E211-BA03-002618B27F8A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0455A25F-0372-E211-9353-003048678F26.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0498C67A-F471-E211-9E4A-0026189438EF.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/061174B4-AB71-E211-ADDA-003048FF86CA.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/066B8A7A-BB71-E211-A2C5-002618FDA259.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/06829A7C-6F71-E211-9921-00304867929E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0682BC2E-FF71-E211-AFCF-0030486790A0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/06BD7444-C471-E211-B6E5-002354EF3BE2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/06FAA92D-6E71-E211-A41E-002618FDA21D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/082F1B4E-7371-E211-82EF-00248C55CC9D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0858206B-D271-E211-89F5-003048678AFA.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/08589611-7071-E211-A4EC-003048679236.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/085DCF4B-F271-E211-B72E-003048FFCBB0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/086318D5-B571-E211-A42D-00261894397D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/087572C7-CF71-E211-994F-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0877391E-F671-E211-B2A2-003048679236.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/08FC68B7-F371-E211-87E1-0026189438A2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0A077994-C771-E211-86E2-003048FFCC18.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0A31BC03-D271-E211-B900-0026189438A5.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0AEA1E9D-F471-E211-AAFD-00261894393D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0C16E3DA-0572-E211-A1FA-003048FFD720.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0C391879-AB71-E211-BF5B-00261894385D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0C6000E6-F271-E211-8C44-0030486792B6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0C8EF937-AC71-E211-9F2D-00304867915A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0C91684D-F271-E211-B366-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0CB5E769-D371-E211-BE21-003048679150.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0CBDE0DC-FF71-E211-8770-002618943956.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0CF0C812-C371-E211-BD9C-00304867D446.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E0A91CA-C771-E211-86C4-0026189437EC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E1A1A8B-0472-E211-A0EF-0026189438FF.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E1FA272-F971-E211-B272-002618943947.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E54F163-F771-E211-8091-003048FFCB8C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E6C45AF-B871-E211-B057-00304867929E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E7B9FB5-D571-E211-B62F-003048679076.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E94E5EB-F471-E211-BC03-002618FDA208.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0E962BB9-B771-E211-92FD-002618943856.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0EDE69C0-D971-E211-B74B-00304867901A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/0EF693B4-0272-E211-8676-00304867908C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/100B5167-F471-E211-8784-00304867918E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1017C895-9871-E211-A094-00304867C1BC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/102F97D3-C371-E211-AFCA-002590596490.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/104D2A14-0472-E211-BF0C-003048678F78.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1050328F-BF71-E211-B054-003048678CA2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/105133E7-A671-E211-9AD9-002618943874.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/106E9789-7471-E211-8382-0026189438CF.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/10AC75AE-6F71-E211-8D6F-003048678F8E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/10CE80D5-BB71-E211-A329-0025905938D4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/120048CB-D471-E211-8E86-00304867D838.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/12094DD7-FB71-E211-97F3-0030486792A8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/122116F1-FB71-E211-9961-002618943918.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/129D8362-0572-E211-8E47-003048FFD770.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/12BD28B4-F371-E211-969E-002354EF3BDD.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/12FBF52C-F271-E211-97E3-002618943924.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/144AF458-0872-E211-96BD-0026189438AB.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1466E243-FF71-E211-9C78-002618943918.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/14703869-0172-E211-AB73-003048679164.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/149DF02A-0072-E211-9FBE-002618943876.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/14B48572-A871-E211-98CB-003048678F74.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/14DF1972-A371-E211-9E32-0030486790A0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/14F8C3B0-F671-E211-9E8A-002618FDA259.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/14FE4B6E-8D71-E211-ACC8-003048679214.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1654CC99-FC71-E211-AE11-00304867903E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/16AAA9CE-FD71-E211-899E-0026189438D8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/16D9C8AD-0272-E211-B6AC-0025905938D4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/16FEB683-6F71-E211-B3D8-0025905822B6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18285F54-0172-E211-9DDC-002618943821.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/184A6CD9-BB71-E211-A272-002618943977.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/186943A3-0A72-E211-988E-0030486791DC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18A194F5-F271-E211-9481-0026189438E1.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18B940C7-C171-E211-BCBF-0026189438FE.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18B94546-A271-E211-A83D-00261894380D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18BD1F2D-0072-E211-A2AD-002618943957.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/18FA56AD-B671-E211-8C14-0025905938A8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A0B74B6-FF71-E211-BCCF-00261894397D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A3EF868-C871-E211-BCAB-003048FFD736.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A40F1B7-0572-E211-8F9E-003048679188.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A4B4C64-DF71-E211-978B-003048678AC0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A7381AC-DA71-E211-999A-0026189438D3.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1A9FF373-C071-E211-9B75-0026189438CE.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1ACB13BE-D971-E211-A2D0-0026189438FD.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1AE36F6A-A571-E211-A9B8-002618943982.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1AF1821D-F471-E211-AE7E-002618943924.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1C116E47-FE71-E211-8F43-003048678BE8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1C17B533-D971-E211-B60D-00261894393F.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1C4064AA-0272-E211-BFC3-002618943800.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1C488108-0272-E211-ABF5-003048679188.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1C7DB095-B671-E211-BEEA-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1CC129EF-F471-E211-B4FA-003048FFD71A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1CD5C5CD-B171-E211-BED3-00304867C0EA.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1CFA21F3-C471-E211-8B28-0025905938B4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1E8FEFCE-F271-E211-8A7C-00304867BF9A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1EE38811-C371-E211-BCF4-003048678B86.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1EED8742-0172-E211-83F6-003048679162.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1EEDA0CA-DD71-E211-9CD7-002590593920.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/1EF5FA85-8E71-E211-B196-003048FFCC2C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/207BCE9B-A471-E211-9748-003048678BAC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/207DA7C8-D371-E211-99DA-0025905938A4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/208042C3-6D71-E211-B2A1-003048678FB8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/20939222-A671-E211-B696-002590593902.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/20F080DE-0772-E211-82F1-002590593872.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/20FB84D0-CA71-E211-920B-003048679164.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/22316FBE-CE71-E211-915B-0025905964A6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/224D27DC-0972-E211-A2D8-00304867918E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/224E3EAD-C171-E211-870F-00304867901A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/22D56B14-FB71-E211-AA4B-003048678E6E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/22DE8CD3-FD71-E211-B171-00248C0BE005.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/22E52F93-F371-E211-B5E3-003048678C06.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/22EDC16B-F471-E211-AAA0-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2414CFFB-BA71-E211-BB20-0026189438CC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/24244947-0372-E211-9A1F-0026189438C0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/24247626-D471-E211-9481-002590593920.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2425266E-B371-E211-B93B-002618FDA28E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/244A8220-D071-E211-9513-002618943886.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2471C7EB-DE71-E211-A38F-0025905938A8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/247BDF7D-F371-E211-B86A-002618943947.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/24A02B31-FC71-E211-9F62-003048678E6E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/24A26DAE-C271-E211-A7E8-003048678B04.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2603A09D-8E71-E211-A556-00259059642A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2606770E-0272-E211-98CB-0026189438ED.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/266D6DDE-F271-E211-AD59-003048678B0E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/26D5F5E4-6F71-E211-B4B0-00304867BFF2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/26E579F5-F271-E211-BE2B-002618943880.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/26EB5865-D571-E211-84E3-002618943876.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28042370-D171-E211-B4A7-002618943984.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2823D69B-BF71-E211-BED2-003048678B38.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28301B5B-F471-E211-8687-003048678FF4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2850B9B1-F371-E211-8E29-002618FDA208.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/285EEA24-B571-E211-B822-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/287639CC-BA71-E211-9300-003048FFCB9E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/287AD495-AC71-E211-A18D-003048678F74.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28831501-B671-E211-8C2B-0026189438FE.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28888CA6-F371-E211-A8F5-003048678ADA.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28B40CCD-B771-E211-9F82-003048FFCC18.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/28E245EE-F271-E211-B581-00304867BED8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2A361880-0172-E211-BD4A-0026189438A5.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2A37B1D0-B771-E211-B026-00261894393C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2A73D9B3-C571-E211-A7D8-002590593920.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2A9C0C67-0172-E211-B485-00261894393C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2AB408A8-AC71-E211-8625-002618943959.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2AE49C39-FC71-E211-AEEA-00261894397D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2AF74A24-0072-E211-AFEF-002618943977.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2C37323B-AC71-E211-B984-0026189438B1.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2C3F36F3-0072-E211-AF65-00304867924E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2C4B0F87-6F71-E211-9FFA-003048678B84.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2C89E1FF-F971-E211-BC16-002618943959.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2CBABC67-F471-E211-BB4F-003048678B36.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2CC0F997-A971-E211-B479-003048678B06.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2CE083F1-B971-E211-A55B-003048FFD732.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2CF22E51-B671-E211-8820-00304867C1B0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2CF25835-C171-E211-8B2F-003048678B86.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E131123-F371-E211-931F-002618FDA262.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E39387E-CD71-E211-9799-00248C0BE014.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E3C3125-FD71-E211-996A-00304867924E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E4ABB8B-F471-E211-B0CC-002354EF3BDD.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E4C5045-F471-E211-90F5-0026189438BD.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E899853-0672-E211-A11A-002618943906.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E90AA3B-0072-E211-906D-003048678E6E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2E940E13-F371-E211-948C-002618943875.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2EC47577-C071-E211-BCD5-003048678FF8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2ED197D0-C471-E211-AB63-0025905964C4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/2EDEBF48-0372-E211-B64F-00304867BEC0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/30322E7A-F471-E211-A842-00261894388D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/30405008-FA71-E211-A744-002590593878.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/306CC9AF-CB71-E211-A820-0026189438B3.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/30CD7666-7071-E211-A867-003048678B38.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/30DED919-0272-E211-AAC4-003048678B0E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/30FF36B3-B871-E211-BDA8-003048679244.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3215DD4E-B971-E211-BB75-002618FDA279.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/324E1CBA-0272-E211-B974-0026189437F0.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/326321AF-0272-E211-98EE-002618943868.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32A8D453-DE71-E211-930C-003048678F78.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32AAAE96-BF71-E211-B113-003048FFCC2C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32ACC7ED-F271-E211-B8BF-002618943937.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32B92D0A-0472-E211-B426-00261894394F.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32E963F0-F371-E211-BC7B-0025905938AA.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/32F65368-7071-E211-BB09-003048678C06.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/340144BA-F371-E211-B94A-003048678B92.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/34015000-CF71-E211-9F55-002618943901.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/34273924-D871-E211-B6C3-002354EF3BE2.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3429A733-F471-E211-9332-002618FDA208.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/344F459F-BF71-E211-816C-003048678F8E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/345D064C-F871-E211-B30D-003048678A80.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/346CC2BE-0572-E211-BE46-002618943807.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/347A8888-F471-E211-9191-003048678B92.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/34F3E1AD-F271-E211-997D-0030486792AC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/36008ECE-F371-E211-BC15-002618943923.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/362B8D85-A471-E211-A52E-002618FDA28E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/362F2C2E-AE71-E211-B616-00304867BFBC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3635154C-C971-E211-9772-00261894387E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3638F463-DE71-E211-953E-0025905964B6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3662524F-FC71-E211-AA29-003048678BE8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/382EEB95-DB71-E211-9B09-00259059649C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/386FFB1D-0C72-E211-BD4E-003048FFCB8C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3881EE43-FF71-E211-9CFC-002618FDA21D.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/388A412D-DD71-E211-966D-003048678BF4.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/38A54F65-F471-E211-AC0E-00248C55CC3C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/38DEA93D-FC71-E211-B40D-003048FFCB8C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/38E52361-F371-E211-AC90-003048FFD71E.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/38E6C0B0-C471-E211-9E1E-003048FFCB6A.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A024D11-A671-E211-BE93-002618943898.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A03E988-CD71-E211-8B7A-0025905938A8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A16311C-0272-E211-A6C6-002618943807.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A2ABC09-F871-E211-883E-003048679244.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A60200E-B471-E211-8030-003048678B94.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A76053E-DE71-E211-A309-00261894393C.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3A7F4F0E-A771-E211-8A3E-0025905964CC.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AA95528-F871-E211-B7E8-00248C55CC62.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AAD6260-B471-E211-9B8B-003048FFCC18.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AB36113-F271-E211-BFD1-003048D15D22.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3ABA601D-0072-E211-8043-0026189438E1.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AC2AC51-AA71-E211-A470-003048678B06.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AD2223E-C471-E211-912F-003048678DD6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3AD5C88F-6F71-E211-94A3-003048679080.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3C2A5333-F471-E211-B0DC-0026189438F8.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3CF1D5AD-C171-E211-A11E-0026189438CB.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3E021448-F471-E211-A563-00304867C034.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3E22C5CD-F271-E211-9D9D-0026189438F6.root',
'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/22Jan2013-v1/20000/3E366582-F371-E211-8932-00261894393F.root',

])


process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

#    fileNames = cms.untracked.vstring(
#'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/200/190/00000/18BB0794-8CDF-E111-B9B0-0025B31E3D3C.root'
 #   ),

    fileNames = filelist,
    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
                            inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )

)


#import EXO.DPAnalysis.skim2012c as fileList
#process.source.fileNames = fileList.fileNames

process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),  
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)   

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.ana = cms.EDAnalyzer('DPAnalysis',
  #  rootFileName     = cms.untracked.string('run2012C.root'),
    rootFileName     = cms.untracked.string('/uscmst1b_scratch/lpc1/3DayLifetime/tnorbert/PhotonData/run2012B1_job.root'),
    triggerName      = cms.vstring('HLT_Photon50_CaloIdVL_IsoL','HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    L1GTSource       = cms.string('L1_SingleEG22'),
    L1Select         = cms.bool( False ),
    isData           = cms.bool( True ),
    cscHaloData      = cms.InputTag("CSCHaloData"),
    staMuons         = cms.InputTag("standAloneMuons"),
    CSCSegmentCollection = cms.InputTag("cscSegments"),
    #DTSegmentCollection = cms.InputTag("dtSegments"),
    DTSegmentCollection = cms.InputTag("dt4DCosmicSegments"),
    muonSource  = cms.InputTag("muonsFromCosmics"),
    trigSource = cms.InputTag("TriggerResults","","HLT"),
    jetSource   = cms.InputTag("ak5PFJets"),
    patJetSource = cms.InputTag("selectedPatJetsPFlow"),
    metSource   = cms.InputTag("pfMet"),
    type1metSource   = cms.InputTag("pfType1CorrectedMet"),
    trackSource = cms.InputTag("generalTracks"),
    electronSource   = cms.InputTag("gsfElectrons"),
    photonSource     = cms.InputTag("myphotons"),
    pvSource         = cms.InputTag("offlinePrimaryVerticesWithBS"),
    beamSpotSource   = cms.InputTag("offlineBeamSpot"),
    EBRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    EERecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    
    BarrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters",""),
    EndcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower",""),
    
    tau                = cms.double( 1000 ), 
    genParticles = cms.InputTag("genParticles"),

    # Set up cuts for physics objects
    # vertex cuts                z   ndof   d0 
    vtxCuts       = cms.vdouble( 99,    0,  99 ),
    # photon cuts                pt   eta  sMajMax,  sMinMin, sMinMax,   dR,  Num  leadingPt  
    photonCuts    = cms.vdouble( 45,  3.5,     99.,      -1.,     99.,   0.0,  1,    45  ),
    # photon isolation           trkR,  ecalSumEt, ecalR, hcalSumEt, hcalR 
    photonIso     = cms.vdouble(  1.,       5.0,   1.,       5.0,   1. ),
    # jet cuts                   pt    eta    dR,  nJets
    jetCuts       = cms.vdouble( 30. , 2.5,  0.3,    0 ),
    metCuts       = cms.vdouble( 0. ),
    # electron cuts              pt  eta  EBIso  EEIso nLostHit  
    electronCuts  = cms.vdouble( 25, 2.5,  0.15,   0.1,      2 ),
    # muon cuts                  pt  eta  Iso  dR   
    muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 ),

)

###########  USE UNCLEANED SUPERCLUSTERS  ######################### 
# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
#process.GlobalTag.globaltag = 'GR_R_53_V18::All'
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag( process.GlobalTag, 'GR_R_53_V10::All' )
process.GlobalTag = GlobalTag( process.GlobalTag, 'GR_R_53_V18::All' )


# to get clustering 
process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder")

process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")
process.uncleanSCRecovered.cleanScCollection=cms.InputTag ("correctedHybridSuperClusters")

# myPhoton sequence
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

import RecoEgamma.EgammaPhotonProducers.photonCore_cfi
import RecoEgamma.EgammaPhotonProducers.photons_cfi

process.myphotonCores=RecoEgamma.EgammaPhotonProducers.photonCore_cfi.photonCore.clone()
process.myphotonCores.scHybridBarrelProducer=cms.InputTag ("uncleanSCRecovered:uncleanHybridSuperClusters")

from RecoEgamma.PhotonIdentification.mipVariable_cfi import *
newMipVariable = mipVariable.clone()
newMipVariable.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
newMipVariable.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import*
newisolationSumsCalculator = isolationSumsCalculator.clone()
newisolationSumsCalculator.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
newisolationSumsCalculator.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

process.myphotons=RecoEgamma.EgammaPhotonProducers.photons_cfi.photons.clone()
process.myphotons.barrelEcalHits=cms.InputTag("reducedEcalRecHitsEB")
process.myphotons.endcapEcalHits=cms.InputTag("reducedEcalRecHitsEE")
process.myphotons.isolationSumsCalculatorSet=newisolationSumsCalculator
process.myphotons.mipVariableSet = newMipVariable
process.myphotons.photonCoreProducer=cms.InputTag("myphotonCores")

process.myPhotonSequence = cms.Sequence(process.myphotonCores+
                                        process.myphotons)
# photonID sequence
from RecoEgamma.PhotonIdentification.photonId_cfi import *
process.myPhotonIDSequence = cms.Sequence(PhotonIDProd)
process.PhotonIDProd.photonProducer=cms.string("myphotons")

process.uncleanPhotons = cms.Sequence(
               process.uncleanSCRecovered *
               process.myPhotonSequence *
               process.myPhotonIDSequence
               )

# PFIso 
#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
#process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

# typeI MET correction 
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

# pat process

# conditions
process.load( "Configuration.Geometry.GeometryIdeal_cff" )
process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.patElectrons.addGenMatch  = False
process.patJets.addGenPartonMatch = False
process.patJets.addGenJetMatch    = False
process.patMETs.addGenMET         = False
process.patMuons.addGenMatch      = False
process.patPhotons.addGenMatch    = False
process.patTaus.addGenMatch       = False
process.patTaus.addGenJetMatch    = False
process.patJetCorrFactors.levels.append( 'L2L3Residual' )

process.out = cms.OutputModule("PoolOutputModule" ,
                fileName = cms.untracked.string( 'patTuple_data.root' ) ,
		outputCommands = cms.untracked.vstring(
			'keep *'
			#               'keep *_cscSegments_*_*'
			#               *patEventContentNoCleaning
			)
																                 )


# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"

usePF2PAT( process
		, runPF2PAT = True
		, jetAlgo   = 'AK5'
		, runOnMC   = False
		, postfix   = postfix
		# for MC
		#, jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
		# for data
		, jetCorrections=('AK5PFchs', ['L2L3Residual'])
	 )



process.p = cms.Path(
                     process.uncleanPhotons*
		     getattr(process,"patPF2PATSequence"+postfix)*
                     process.producePFMETCorrections *
                     process.ana
                    )

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False

# enable delta beta correction for muon selection in PF2PAT?
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

process.out.outputCommands.extend( [ 'drop *_*_*_*' ] )

