<?xml version="1.0" encoding="UTF-8"?>
<cs:chemicalSemantics xmlns:cs="http://chemicalsemantics.com/csx" version="1.0">
    <cs:molecularPublication xmlns:cs="http://chemicalsemantics.com/csx" xmlns:dcterms="http://purl.org/dc/terms/" >
        <dcterms:title>default title</dcterms:title>
        <dcterms:abstract>default abstract</dcterms:abstract>
        <dcterms:publisher>default publisher</dcterms:publisher>
        <cs:author type="cs:corresponding">
            <dcterms:creator>Wang</dcterms:creator>
            <cs:organization>default organization</cs:organization>
            <cs:email>name@organization.com</cs:email>
        </cs:author>
        <cs:sourcePackage>
            <cs:name>DALTON</cs:name>
            <cs:version>2015.0</cs:version>
        </cs:sourcePackage>
        <cs:tags>DALTON</cs:tags>
        <cs:status>default status</cs:status>
        <cs:visibility></cs:visibility>
        <cs:category>1</cs:category>
        <cs:key>1</cs:key>
    </cs:molecularPublication>
    <cs:molecularSystem>
        <cs:systemTemperature unit="cs:kelvin">0.0</cs:systemTemperature>
        <cs:systemCharge>0</cs:systemCharge>
        <cs:systemMultiplicity>1</cs:systemMultiplicity>
        <cs:molecule atomCount="3" id="m1">
            <cs:atom id="a1">
                <cs:elementSymbol>O</cs:elementSymbol>
                <cs:atomMass>15.994914619559999</cs:atomMass>
                <cs:formalAtomCharge>0</cs:formalAtomCharge>
                <cs:calculatedAtomCharge>0.</cs:calculatedAtomCharge>
                <cs:xCoord3D unit="cs:angstrom">0.0</cs:xCoord3D>
                <cs:yCoord3D unit="cs:angstrom">0.0</cs:yCoord3D>
                <cs:zCoord3D unit="cs:angstrom">-0.0352847554102</cs:zCoord3D>
                <cs:basisSet>cs:STO-3G</cs:basisSet>
            </cs:atom>
            <cs:atom id="a2">
                <cs:elementSymbol>H</cs:elementSymbol>
                <cs:atomMass>1.00782503207</cs:atomMass>
                <cs:formalAtomCharge>0</cs:formalAtomCharge>
                <cs:calculatedAtomCharge>0.</cs:calculatedAtomCharge>
                <cs:xCoord3D unit="cs:angstrom">0.0</cs:xCoord3D>
                <cs:yCoord3D unit="cs:angstrom">0.418393547929</cs:yCoord3D>
                <cs:zCoord3D unit="cs:angstrom">0.279997320777</cs:zCoord3D>
                <cs:basisSet>cs:STO-3G</cs:basisSet>
            </cs:atom>
            <cs:atom id="a3">
                <cs:elementSymbol>H</cs:elementSymbol>
                <cs:atomMass>1.00782503207</cs:atomMass>
                <cs:formalAtomCharge>0</cs:formalAtomCharge>
                <cs:calculatedAtomCharge>0.</cs:calculatedAtomCharge>
                <cs:xCoord3D unit="cs:angstrom">0.0</cs:xCoord3D>
                <cs:yCoord3D unit="cs:angstrom">-0.418393547929</cs:yCoord3D>
                <cs:zCoord3D unit="cs:angstrom">0.279997320777</cs:zCoord3D>
                <cs:basisSet>cs:STO-3G</cs:basisSet>
            </cs:atom>
        </cs:molecule>
    </cs:molecularSystem>
    <cs:molecularCalculation>
        <cs:quantumMechanics>
            <cs:singleReferenceState>
                <cs:singleDeterminant>
                    <cs:mp2 spinType="cs:RHF" basisSet="bse:STO-3G" methodology="cs:normal">
                        <cs:energies unit="cs:eV">
                            <cs:energy type="cs:totalPotential">-1998.57754348</cs:energy>
                            <cs:energy type="cs:correlation">-0.320694682331</cs:energy>
                        </cs:energies>
                    </cs:mp2>
                </cs:singleDeterminant>
            </cs:singleReferenceState>
        </cs:quantumMechanics>
    </cs:molecularCalculation>
</cs:chemicalSemantics>
