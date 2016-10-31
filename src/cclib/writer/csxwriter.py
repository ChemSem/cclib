#the generation of a CSX file from a quantum chemistry package output file
#using cclib as a parser

from . import filewriter

class CSX(filewriter.Writer):

    def __init__(self, ccdata, splitfiles=False,
            firstgeom=False, lastgeom=True, allgeom=False,
            *args, **kwargs):
        
    # Call the __init__ method of the superclass
        super(CSX, self).__init__(ccdata, *args, **kwargs)

        self.do_firstgeom = firstgeom
        self.do_lastgeom = lastgeom
        self.do_allgeom = allgeom

        self.generate_csx()

    def generate_csx(self):
        import os
        import sys
        import math
        #import glob
        #import openbabel
        import chemElement
        import csx2_api as api
        #from cclib.parser import ccopen
        #from cclib.writer import ccwrite


        #inputfile = sys.argv[1]
        #myFile = ccopen(inputfile)
        fileName, fileExt = os.path.splitext(self.jobfilename)
        #fileName = str(self.jobfilename)
        csxfile = open( fileName + '.csx', 'w')
        data = self.ccdata
        #cmlfile = open( fileName + '.cml', 'w')
        #xyz = ccwrite(data, 'cml', cmlfile)
        atomNum = data.natom
        if (hasattr(data, 'basisname')):
            basisName = data.basisname
        else:
            basisName = 'none'
        molCharge = data.charge
        molMulti = data.mult
        wfnRestricted = False
        hasOrb = True if (hasattr(data, 'mocoeffs')) else False
        hasOrbSym = True if (hasattr(data, 'mosyms')) else False
        hasFreq = True if (hasattr(data, 'vibfreqs')) else False
        hasProp = True if (hasattr(data, 'moments')) else False
        hasPolar = True if (hasattr(data, 'polar')) else False
        hasNMR = True if (hasattr(data, 'nmriso')) else False
        hasElec = True if (hasattr(data, 'etoscs')) else False
        hasCfunctional = True if (hasattr(data, 'cfunctional')) else False
        hasCollection = True if (hasattr(data, 'ircpnt')) else False
        molSpin = data.spintype if (hasattr(data, 'spintype')) else 'RHF'
        if molMulti == 1 :
            wfnRestricted = True
        else:
            if molSpin == 'ROHF':
                wfnRestricted = True
            else:
                molSpin = 'UHF'
        calcType = data.theory
        molEE = data.scfenergies[-1]
        #Wavefunction
        if hasOrb:
            if wfnRestricted :
                orbEne = data.moenergies
                orbEString = ' '.join(str(x) for x in orbEne[0])
                orbNum = data.nmo
                orbOcc = []
                if hasOrbSym:
                    orbSym = data.mosyms
                    for x in orbSym[0]:
                        if '\"' in x:
                            xindex = orbSym[0].index(x)
                            x = x[:-1]+'\''+'\''
                            orbSym[0][xindex] = x
                    orbSymString = ' '.join( x for x in orbSym[0])
                for iorb in range (orbNum):
                    if molSpin == "ROHF":
                        if iorb < int(data.homos[0]):
                            elecNum = 2
                        elif iorb == int(data.homos[0]):
                            elecNum = 1
                        else:
                            elecNum = 0
                    else:
                        elecNum = 0 if iorb > int(data.homos) else 2
                    orbOcc.append(elecNum)
                orbOccString = ' '.join(str(x) for x in sorted(orbOcc,reverse=True))
                wfn1 = api.waveFunctionType(orbitalCount=orbNum, \
                        orbitalOccupancies=orbOccString)
                if hasOrbSym:
                    wfn1.set_orbitalSymmetry(orbSymString)
                orbe1 = api.stringArrayType(unit='u:ElectronVolt')
                orbe1.set_valueOf_(orbEString)
                orbs1 = api.orbitalsType()
                for orbArray in data.mocoeffs:
                    for iorb in range(orbNum):
                        orbCaString = ' '.join(str(x) for x in orbArray[iorb])
                        orb1 = api.stringArrayType(id=iorb+1)
                        orb1.set_valueOf_(orbCaString)
                        orbs1.add_orbital(orb1)
                wfn1.set_orbitals(orbs1)
                wfn1.set_orbitalEnergies(orbe1)
            else:
                orbNum = data.nmo
                orbCaEne = data.moenergies[0][:]
                orbCaEString = ' '.join(str(x) for x in orbCaEne)
                orbCbEne = data.moenergies[1][:]
                orbCbEString = ' '.join(str(x) for x in orbCbEne)
                orbCaOcc = []
                orbCbOcc = []
                for iorb in range (orbNum):
                    elecCa = 1 if orbCaEne[iorb] < 0.0 else 0
                    orbCaOcc.append(elecCa)
                    elecCb = 1 if orbCbEne[iorb] < 0.0 else 0
                    orbCbOcc.append(elecCb)
                orbCaOccString = ' '.join(str(x) for x in sorted(orbCaOcc,reverse=True))
                orbCbOccString = ' '.join(str(x) for x in sorted(orbCbOcc,reverse=True))
                orbCaSym = data.mosyms[0][:]
                orbCaSymString = ' '.join( x for x in orbCaSym)
                orbCbSym = data.mosyms[1][:]
                orbCbSymString = ' '.join( x for x in orbCbSym)
                wfn1 = api.waveFunctionType(orbitalCount=orbNum)
                orbe1 = api.stringArrayType(unit='u:ElectronVolt')
                orbe1.set_valueOf_(orbCaEString)
                orbs1 = api.orbitalsType()
                alphaOrb = data.mocoeffs[0][:]
                for iorb in range(orbNum):
                    orbCaString = ' '.join(str(x) for x in alphaOrb[iorb])
                    orb1 = api.stringArrayType(id=iorb+1)
                    orb1.set_valueOf_(orbCaString)
                    orbs1.add_orbital(orb1)
                wfn1.set_alphaOrbitals(orbs1)
                wfn1.set_alphaOrbitalEnergies(orbe1)
                wfn1.set_alphaOrbitalOccupancies(orbCaOccString)
                wfn1.set_alphaOrbitalSymmetry(orbCaSymString)
                orbe2 = api.stringArrayType(unit='u:ElectronVolt')
                orbe2.set_valueOf_(orbCbEString)
                orbs2 = api.orbitalsType()
                betaOrb = data.mocoeffs[1][:]
                for iorb in range(orbNum):
                    orbCbString = ' '.join(str(x) for x in betaOrb[iorb])
                    orb2 = api.stringArrayType(id=iorb+1)
                    orb2.set_valueOf_(orbCbString)
                    orbs2.add_orbital(orb2)
                wfn1.set_betaOrbitals(orbs2)
                wfn1.set_betaOrbitalEnergies(orbe2)
                wfn1.set_betaOrbitalOccupancies(orbCbOccString)
                wfn1.set_betaOrbitalSymmetry(orbCbSymString)

        #vibrational frequency
        if hasFreq:
            molFreqNum = len(data.vibfreqs)
            frqString = ' '.join(str(x) for x in data.vibfreqs)
            if hasattr(data, 'vibirs'):
                intString = ' '.join(str(x) for x in data.vibirs)
            else:
                intString = ' '.join(str(x) for x in [1.0]*molFreqNum)
            vib1 = api.vibAnalysisType(vibrationCount=molFreqNum)
            freq1 = api.stringArrayType(unit="gc:RecipricalCentimeter")
            freq1.set_valueOf_(frqString)
            vib1.set_frequencies(freq1)
            irint1 = api.stringArrayType()
            irint1.set_valueOf_(intString)
            vib1.set_irIntensities(irint1)
            norms1 = api.normalModesType()
            normMdString = []
            for ifrq in range(molFreqNum):
                normM = []
                for iatm in range(atomNum):
                    for ixyz in range(3):
                        normM.append(data.vibdisps[ifrq][iatm][ixyz])
                normMdString.append(' '.join(str(x) for x in normM))
                norm1 = api.normalModeType(id=ifrq+1)
                norm1.set_valueOf_(normMdString[ifrq])
                norms1.add_normalMode(norm1)
            vib1.set_normalModes(norms1)

        if hasProp or hasPolar:
            prop1 = api.propertiesType()
        #dipole moments information
        if hasProp:
        #    au2db = 2.541766
            if len(data.moments) > 1:
                molDipoleX = data.moments[1][0]
                molDipoleY = data.moments[1][1]
                molDipoleZ = data.moments[1][2]
                molDipoleTot = math.sqrt(molDipoleX*molDipoleX+molDipoleY*molDipoleY+molDipoleZ*molDipoleZ)
                sprop1 = api.propertyType(name='dipoleMomentX',unit='u:Debye',moleculeId='m1')
                sprop1.set_valueOf_(molDipoleX)
                sprop2 = api.propertyType(name='dipoleMomentY',unit='u:Debye',moleculeId='m1')
                sprop2.set_valueOf_(molDipoleY)
                sprop3 = api.propertyType(name='dipoleMomentZ',unit='u:Debye',moleculeId='m1')
                sprop3.set_valueOf_(molDipoleZ)
                sprop4 = api.propertyType(name='dipoleMomentAverage',unit='u:Debye',moleculeId='m1')
                sprop4.set_valueOf_(molDipoleTot)
                prop1.add_systemProperty(sprop1)
                prop1.add_systemProperty(sprop2)
                prop1.add_systemProperty(sprop3)
                prop1.add_systemProperty(sprop4)
        #polarizability information
        if hasPolar:
            polarXX = data.polar[0][0]
            polarYY = data.polar[0][1]
            polarZZ = data.polar[0][2]
            polarAvg = (polarXX+polarYY+polarZZ)/3.0
            sprop5 = api.propertyType(name='polarizabilityXX',unit='u:Angstrom3',moleculeId='m1')
            sprop5.set_valueOf_(polarXX)
            sprop6 = api.propertyType(name='polarizabilityYY',unit='u:Angstrom3',moleculeId='m1')
            sprop6.set_valueOf_(polarYY)
            sprop7 = api.propertyType(name='polarizabilityZZ',unit='u:Angstrom3',moleculeId='m1')
            sprop7.set_valueOf_(polarZZ)
            sprop8 = api.propertyType(name='polarizabilityAverage',unit='u:Angstrom3')
            sprop8.set_valueOf_(polarAvg)
            prop1.add_systemProperty(sprop5)
            prop1.add_systemProperty(sprop6)
            prop1.add_systemProperty(sprop7)
            prop1.add_systemProperty(sprop8)
        #NMR chemical shielding information
        if hasNMR:
            prop2 = api.propertiesType()
            for iatm in range(atomNum):
                aprop1 = api.propertyType(atomId='a'+str(iatm+1), moleculeId='m1', \
                        name='nmrShieldingIsotropic',unit='gc:PartsPerMillion')
                aprop1.set_valueOf_(data.nmriso[iatm])
                aprop2 = api.propertyType(atomId='a'+str(iatm+1), moleculeId='m1', \
                        name='nmrShieldingAnisotropy',unit='gc:PartsPerMillion')
                aprop2.set_valueOf_(data.nmranis[iatm])
                prop2.add_atomProperty(aprop1)
                prop2.add_atomProperty(aprop2)

        #Electronic transition information
        if hasElec:
            transStr = ' '.join(str(x) for x in data.etenergies)
            oscilStr = ' '.join(str(x) for x in data.etoscs)
            elec1 = api.elecSpectraType(transitionCount=len(data.etenergies))
            trans1 = api.stringArrayType(unit="gc:RecipricalCentimeter")
            trans1.set_valueOf_(transStr)
            oscil1 = api.stringArrayType()
            oscil1.set_valueOf_(oscilStr)
            elec1.set_electronicTransitions(trans1)
            elec1.set_oscillatorStrength(oscil1)

        #Start to generate CSX elements
        cs1 = api.csType(version='2.0')

        #molecular publication section
        mp1 = api.mpubType(title='default title', \
                abstract='default abstract', \
                publisher='default publisher', \
                status='default status', \
                category=1, \
                visibility=0, \
                tags=data.package, \
                key=1 )
        source1 = api.sourcePackageType(name=data.package, version=data.version)
        mp1.set_sourcePackage(source1)
        ath1 = api.authorType(creator='Wang', \
                type_='gc:CorrespondingAuthor', \
                organization='default organization', \
                email='name@organization.com')
        mp1.add_author(ath1)
        cs1.set_molecularPublication(mp1)

        if hasCollection:
            mc1 = api.mcolType(type_='IRC', description='Internal Reaction Coordinates')
            ircpnt = int(data.ircpnt)
            if ( ircpnt != len(data.ircenergies) ):
                ircpnt = len(data.ircenergies)-1
            for ipnt in range(ircpnt):
                entry1 = api.entryType(id='pnt'+str(ipnt+1))
                param1 = api.propertyType(id='pnt'+str(ipnt+1)+'prm1', name='reaction coordinates')
                param1.set_valueOf_(data.irccoords[ipnt])
                entry1.add_parameter(param1)
                sys1 = api.itemType(ref='s'+str(ipnt+1))
                calc1 = api.itemType(ref='c'+str(ipnt+1))
                res1 = api.resType(id='pnt'+str(ipnt+1)+'r1', \
                        name='total electronic energy for point '+str(ipnt+1))
                val1 = api.valType(ref='e'+str(ipnt+1)+'_'+'total_energy')
                res1.add_value(val1)
                entry1.add_result(res1)
                entry1.set_system(sys1)
                entry1.set_calculation(calc1)
                mc1.add_entry(entry1)
            cs1.add_molecularCollection(mc1)

            for ipnt in range(ircpnt+1):
                #molecular system section
                msys = api.msysType(systemCharge=molCharge, \
                        systemMultiplicity=molMulti, id='s'+str(ipnt))
                temp1 = api.dataWithUnitsType(unit='u:Kelvin')
                temp1.set_valueOf_(0.0)
                msys.set_systemTemperature(temp1)
                if ipnt == 0:
                    mol = api.moleculeType(id='m1',atomCount=atomNum)
                else:
                    mol = api.moleculeType(ref='m1',atomCount=atomNum)
                if hasattr(data, "atomcharges"):
                    atmCharge = data.atomcharges["mulliken"]
                else:
                    atmCharge = [0]*atomNum
                for iatm in range(atomNum):
                    #   xCoord = float(data.atomcoords[iatm,0])
                    xCoord = data.atomcoords[ipnt,iatm,0]
                    yCoord = data.atomcoords[ipnt,iatm,1]
                    zCoord = data.atomcoords[ipnt,iatm,2]
                    xCoord1 = api.dataWithUnitsType(unit='u:Angstrom')
                    xCoord1.set_valueOf_(xCoord)
                    yCoord1 = api.dataWithUnitsType(unit='u:Angstrom')
                    yCoord1.set_valueOf_(yCoord)
                    zCoord1 = api.dataWithUnitsType(unit='u:Angstrom')
                    zCoord1.set_valueOf_(zCoord)
                    atomicNum = data.atomnos[iatm]
                    if ipnt == 0:
                        atm = api.atomType(id='a'+str(iatm+1), elementSymbol=chemElement.z2elm[atomicNum], \
                                atomMass=chemElement.z2mass[atomicNum], \
                                xCoord3D=xCoord1, \
                                yCoord3D=yCoord1, \
                                zCoord3D=zCoord1, \
                                basisSet='bse:'+basisName, \
                                calculatedAtomCharge=atmCharge[iatm], \
                                formalAtomCharge=0)
                    else:
                        atm = api.atomType(ref='a'+str(iatm+1), elementSymbol=chemElement.z2elm[atomicNum], \
                                atomMass=chemElement.z2mass[atomicNum], \
                                xCoord3D=xCoord1, \
                                yCoord3D=yCoord1, \
                                zCoord3D=zCoord1, \
                                basisSet='bse:'+basisName, \
                                calculatedAtomCharge=atmCharge[iatm], \
                                formalAtomCharge=0)
                    mol.add_atom(atm)
                msys.add_molecule(mol)
                cs1.add_molecularSystem(msys)

            for ipnt in range(ircpnt):
                #molCalculation section
                sd_wfn_method = ['HF', 'DFT', 'MP2', 'MP3', 'MP4', 'AM1', 'PM3', 'PM6']
                md_wfn_method = ['CCD', 'CCSD', 'CCSD-T', 'CIS', 'CISD', 'FCI', 'QCISD', 'QCISD-T']
                mr_wfn_method = ['CASSCF', 'CASPT2', 'RASSCF', 'RASPT2', 'GVB', 'MCSCF', 'MRCC', 'MRCI']
                mc1 = api.mcalType(id='c'+str(ipnt+1), inputsystem='s'+str(ipnt), \
                        outputsystem='s'+str(ipnt+1))
                qm1 = api.qmCalcType()
                if calcType in mr_wfn_method:
                    mrs1 = api.mrsMethodType()
                    mrsmd1 = api.mrsmdMethodType()
                    if calcType == 'GVB':
                        gvb1 = api.resultType(methodology='gc:normal', spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName, pairCount='2')
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(molEE))
                        ene1.add_energy(ee_ene1)
                        gvb1.set_energies(ene1)
                        mrsmd1.set_gvb(gvb1)
                    elif calcType == 'CASSCF':
                        casscf1  = api.resultType(methodology='gc:normal', spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(molEE))
                        ene1.add_energy(ee_ene1)
                        casscf1.set_energies(ene1)
                        mrsmd1.set_casscf(casscf1)
                    mrs1.set_multipleDeterminant(mrsmd1)
                    qm1.set_multipleReferenceState(mrs1)
                else:
                    srs1 = api.srsMethodType()
                    if calcType in sd_wfn_method:
                        sdm1 = api.srssdMethodType()
                        #SCF
                        if (calcType == 'HF'):
                            scf1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                    basisSet='bse:'+basisName)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(id='e'+str(ipnt+1)+'_'+'total_energy', \
                                    type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(data.ircenergies[ipnt])
                            ene1.add_energy(ee_ene1)
                            scf1.set_energies(ene1)
                            sdm1.set_abinitioScf(scf1)
                        #DFT
                        elif (calcType == 'DFT'):
                            if hasElec:
                                dft1 = api.resultType(methodology='gc:tddft',spinType='gc:'+molSpin, \
                                        basisSet='bse:'+basisName, dftFunctional='gc:'+data.functional)
                            else:
                                dft1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                        basisSet='bse:'+basisName, dftFunctional='gc:'+data.functional)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(id='e'+str(ipnt+1)+'_'+'total_energy', \
                                    type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(data.ircenergies[ipnt])
                            ene1.add_energy(ee_ene1)
                            dft1.set_energies(ene1)
                            sdm1.set_dft(dft1)
                        #MP2
                        elif (calcType == 'MP2'):
                            mp21 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                    basisSet='bse:'+basisName)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(id='e'+str(ipnt+1)+'_'+'total_energy', \
                                    type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(data.ircenergies[ipnt])
                            ene1.add_energy(ee_ene1)
                            mp21.set_energies(ene1)
                            sdm1.set_mp2(mp21)
                        #Semiempirical methods
                        elif (calcType == 'AM1' or calcType == 'PM3' or calcType == 'PM6'):
                            sem1 = api.resultType(methodology=calcType,spinType='gc:'+molSpin)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(float(molEE))
                            hof_ene1 = api.energyType(type_='gc:heatofformation')
                            hof_ene1.set_valueOf_(float(data.hofenergies[-1]))
                            ene1.add_energy(ee_ene1)
                            ene1.add_energy(hof_ene1)
                            sem1.set_energies(ene1)
                            sdm1.set_semiEmpiricalScf(sem1)
                        else:
                            print ('The current CSX does not support this method')

                        srs1.set_singleDeterminant(sdm1)

                    if calcType in md_wfn_method:
                        mdm1 = api.srsmdMethodType()
                        if (calcType == 'CCSD'):
                            ccsd1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                    basisSet='bse:'+basisName)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(float(data.ccenergies[0]))
                            ce_ene1 = api.energyType(type_='gc:correlation')
                            ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                            ene1.add_energy(ee_ene1)
                            ene1.add_energy(ce_ene1)
                            ccsd1.set_energies(ene1)
                            mdm1.set_ccsd(ccsd1)
                        elif (calcType == 'CCSD-T'):
                            ccsd_t1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                    basisSet='bse:'+basisName)
                            ene1 = api.energiesType(unit='u:ElectronVolt')
                            ee_ene1 = api.energyType(type_='gc:totalPotential')
                            ee_ene1.set_valueOf_(float(data.ccenergies[-1]))
                            ce_ene1 = api.energyType(type_='gc:correlation')
                            ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                            ene1.add_energy(ee_ene1)
                            ene1.add_energy(ce_ene1)
                            ccsd_t1.set_energies(ene1)
                            mdm1.set_ccsd-t(ccsd_t1)
                        else:
                            print ('The current CSX does not support this method')
                        srs1.set_multipleDeterminant(mdm1)
                    qm1.set_singleReferenceState(srs1)
                mc1.set_quantumMechanics(qm1)
                cs1.add_molecularCalculation(mc1)
        else:
            #molecular system section
            ms1 = api.msysType(systemCharge=molCharge, \
                   systemMultiplicity=molMulti, id='s1')
            temp1 = api.dataWithUnitsType(unit='u:Kelvin')
            temp1.set_valueOf_(0.0)
            ms1.set_systemTemperature(temp1)
            mol1 = api.moleculeType(id='m1',atomCount=atomNum)
            if hasattr(data, "atomcharges"):
                atmCharge = data.atomcharges["mulliken"]
            else:
                atmCharge = [0]*atomNum
            #obmol1 = openbabel.OBMol()
            for iatm in range(atomNum):
                #   xCoord = float(data.atomcoords[iatm,0])
                xCoord = data.atomcoords[-1,iatm,0]
                yCoord = data.atomcoords[-1,iatm,1]
                zCoord = data.atomcoords[-1,iatm,2]
                xCoord1 = api.floatWithUnitType(unit='u:Angstrom')
                xCoord1.set_valueOf_(xCoord)
                yCoord1 = api.floatWithUnitType(unit='u:Angstrom')
                yCoord1.set_valueOf_(yCoord)
                zCoord1 = api.floatWithUnitType(unit='u:Angstrom')
                zCoord1.set_valueOf_(zCoord)
                atomicNum = data.atomnos[iatm]
                atm = api.atomType(id='a'+str(iatm+1), elementSymbol=chemElement.z2elm[atomicNum], \
                        atomMass=chemElement.z2mass[atomicNum], \
                        xCoord3D=xCoord1, \
                        yCoord3D=yCoord1, \
                        zCoord3D=zCoord1, \
                        basisSet='bse:'+basisName, \
                        calculatedAtomCharge=atmCharge[iatm], \
                        formalAtomCharge=0)
                mol1.add_atom(atm)
            ms1.add_molecule(mol1)
            cs1.add_molecularSystem(ms1)

            #molCalculation section
            sd_wfn_method = ['HF', 'DFT', 'MP2', 'MP3', 'MP4', 'AM1', 'PM3', 'PM6']
            md_wfn_method = ['CCD', 'CCSD', 'CCSD-T', 'CIS', 'CISD', 'FCI', 'QCISD', 'QCISD-T']
            mr_wfn_method = ['CASSCF', 'CASPT2', 'RASSCF', 'RASPT2', 'GVB', 'MCSCF', 'MRCC', 'MRCI']
            mc1 = api.mcalType(id='c1', inputsystem='s1')
            qm1 = api.qmCalcType()
            if calcType in mr_wfn_method:
                mrs1 = api.mrsMethodType()
                mrsmd1 = api.mrsmdMethodType()
                if calcType == 'GVB':
                    gvb1 = api.resultType(methodology='gc:normal', spinType='gc:'+molSpin, \
                            basisSet='bse:'+basisName, pairCount='2')
                    ene1 = api.energiesType(unit='u:ElectronVolt')
                    ee_ene1 = api.energyType(type_='gc:totalPotential')
                    ee_ene1.set_valueOf_(float(molEE))
                    ene1.add_energy(ee_ene1)
                    gvb1.set_energies(ene1)
                    mrsmd1.set_gvb(gvb1)
                elif calcType == 'CASSCF':
                    casscf1  = api.resultType(methodology='gc:normal', spinType='gc:'+molSpin, \
                            basisSet='bse:'+basisName)
                    ene1 = api.energiesType(unit='u:ElectronVolt')
                    ee_ene1 = api.energyType(type_='gc:totalPotential')
                    ee_ene1.set_valueOf_(float(molEE))
                    ene1.add_energy(ee_ene1)
                    casscf1.set_energies(ene1)
                    mrsmd1.set_casscf(casscf1)
                elif calcType == 'MCSCF':
                    mcscf1  = api.resultType(methodology='gc:normal', spinType='gc:'+molSpin, \
                            basisSet='bse:'+basisName)
                    ene1 = api.energiesType(unit='u:ElectronVolt')
                    ee_ene1 = api.energyType(type_='gc:totalPotential')
                    ee_ene1.set_valueOf_(float(molEE))
                    ene1.add_energy(ee_ene1)
                    mcscf1.set_energies(ene1)
                    mrsmd1.set_mcscf(mcscf1)
                mrs1.set_multipleDeterminant(mrsmd1)
                qm1.set_multipleReferenceState(mrs1)
            else:
                srs1 = api.srsMethodType()
                if calcType in sd_wfn_method:
                    sdm1 = api.srssdMethodType()
                    #SCF
                    if (calcType == 'HF'):
                        scf1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(molEE))
                        ene1.add_energy(ee_ene1)
                        scf1.set_energies(ene1)
                        if hasOrb:
                            scf1.set_waveFunction(wfn1)
                        if hasProp or hasPolar:
                            scf1.set_properties(prop1)
                        if hasNMR:
                            scf1.set_properties(prop2)
                        if hasFreq:
                            scf1.set_vibrationalAnalysis(vib1)
                        if hasElec:
                            scf1.set_electronicSpectra(elec1)
                        sdm1.set_abinitioScf(scf1)
                    #DFT
                    elif (calcType == 'DFT'):
                        if hasElec:
                            dft1 = api.resultType(methodology='gc:tddft',spinType='gc:'+molSpin, \
                                    basisSet='bse:'+basisName, dftFunctional='gc:'+data.functional)
                        else:
                            if hasCfunctional:
                                dft1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                        basisSet='bse:'+basisName, \
                                        exchangeFunctional='gc:'+data.xfunctional, \
                                        correlationFunctional='gc:'+data.cfunctional)
                            else:
                                dft1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                        basisSet='bse:'+basisName, dftFunctional='gc:'+data.functional)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(molEE))
                        ene1.add_energy(ee_ene1)
                        dft1.set_energies(ene1)
                        if hasOrb:
                            dft1.set_waveFunction(wfn1)
                        if hasProp:
                            dft1.set_properties(prop1)
                        if hasNMR:
                            dft1.set_properties(prop2)
                        if hasFreq:
                            dft1.set_vibrationalAnalysis(vib1)
                        if hasElec:
                            dft1.set_electronicSpectra(elec1)
                        sdm1.set_dft(dft1)
                    #MP2
                    elif (calcType == 'MP2'):
                        mp21 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(data.mpenergies[-1]))
                        ce_ene1 = api.energyType(type_='gc:correlation')
                        ce_ene1.set_valueOf_(float(data.mpenergies[-1])-float(molEE))
                        ene1.add_energy(ee_ene1)
                        ene1.add_energy(ce_ene1)
                        mp21.set_energies(ene1)
                        if hasOrb:
                            mp21.set_waveFunction(wfn1)
                        if hasProp:
                            mp21.set_properties(prop1)
                        if hasNMR:
                            mp21.set_properties(prop2)
                        if hasFreq:
                            mp21.set_vibrationalAnalysis(vib1)
                        sdm1.set_mp2(mp21)
                    #Semiempirical methods
                    elif (calcType == 'AM1' or calcType == 'PM3' or calcType == 'PM6'):
                        sem1 = api.resultType(methodology=calcType,spinType='gc:'+molSpin)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(molEE))
                        hof_ene1 = api.energyType(type_='gc:heatofformation')
                        hof_ene1.set_valueOf_(float(data.hofenergies[-1]))
                        ene1.add_energy(ee_ene1)
                        ene1.add_energy(hof_ene1)
                        sem1.set_energies(ene1)
                        if hasOrb:
                            sem1.set_waveFunction(wfn1)
                        if hasProp:
                            sem1.set_properties(prop1)
                        if hasNMR:
                            sem1.set_properties(prop2)
                        if hasFreq:
                            sem1.set_vibrationalAnalysis(vib1)
                        sdm1.set_semiEmpiricalScf(sem1)
                    else:
                        print ('The current CSX does not support this method')

                    srs1.set_singleDeterminant(sdm1)

                if calcType in md_wfn_method:
                    mdm1 = api.srsmdMethodType()
                    if (calcType == 'CCSD'):
                        ccsd1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(data.ccenergies[0]))
                        ce_ene1 = api.energyType(type_='gc:correlation')
                        ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                        ene1.add_energy(ee_ene1)
                        ene1.add_energy(ce_ene1)
                        ccsd1.set_energies(ene1)
                        if hasOrb:
                            ccsd1.set_waveFunction(wfn1)
                        if hasProp:
                            ccsd1.set_properties(prop1)
                        if hasNMR:
                            ccsd1.set_properties(prop2)
                        if hasFreq:
                            ccsd1.set_vibrationalAnalysis(vib1)
                        mdm1.set_ccsd(ccsd1)
                    elif (calcType == 'CCSD-T'):
                        ccsd_t1 = api.resultType(methodology='gc:normal',spinType='gc:'+molSpin, \
                                basisSet='bse:'+basisName)
                        ene1 = api.energiesType(unit='u:ElectronVolt')
                        ee_ene1 = api.energyType(type_='gc:totalPotential')
                        ee_ene1.set_valueOf_(float(data.ccenergies[-1]))
                        ce_ene1 = api.energyType(type_='gc:correlation')
                        ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                        ene1.add_energy(ee_ene1)
                        ene1.add_energy(ce_ene1)
                        ccsd_t1.set_energies(ene1)
                        if hasOrb:
                            ccsd_t1.set_waveFunction(wfn1)
                        if hasProp:
                            ccsd_t1.set_properties(prop1)
                        if hasNMR:
                            ccsd_t1.set_properties(prop2)
                        if hasFreq:
                            ccsd_t1.set_vibrationalAnalysis(vib1)
                        mdm1.set_ccsd-t(ccsd_t1)
                    else:
                        print ('The current CSX does not support this method')
                    srs1.set_multipleDeterminant(mdm1)
                qm1.set_singleReferenceState(srs1)
            mc1.set_quantumMechanics(qm1)
            cs1.add_molecularCalculation(mc1)

        csxfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        cs1.export(csxfile, 0)
        csxfile.close()


