/*

AIM:
~~~

Plot the distribution of angles between various vectors.

The vectors are defined as follows,
(i) From head group N to tail end C
(ii) From C(n) to C(n+x), where x is the variable less than tail length
(iii) Between every bonds within a cutoff distance
(iv) Normal vector to Au surface and head N to tail C
*/

/*
INPUTS:
~~~~~~

(i) input dump file
(ii) input data file
(iii) 

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>

#define PI 3.14159
#define ORDERPARAMETERBINDISTANCE 1.0

typedef struct trajectory
{
	int atomID, atomType, molType, ix, iy, iz;
	float x, y, z;
	int isEndGroup;
} TRAJECTORY;

typedef struct vector
{
	float x1, y1, z1;
	float x2, y2, z2;
	float xc, yc, zc;
} VECTOR;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct rdf
{
	float rlo, rhi, gofr;
} RDF;

typedef struct stats
{
	float average, standardDeviation;
} STATS;

typedef struct orderParameterBins
{
	long double orderParameter, rlo, rhi, count;
} ORDERPARAMETER_BINS;

SIMULATION_BOUNDARY readDumpBoundary (FILE *file_dump, SIMULATION_BOUNDARY boundary)
{
	rewind (file_dump);
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}

	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.xlo, &boundary.xhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.ylo, &boundary.yhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.zlo, &boundary.zhi);
	rewind (file_dump);

	boundary.xLength = boundary.xhi - boundary.xlo;
	boundary.yLength = boundary.yhi - boundary.ylo;
	boundary.zLength = boundary.zhi - boundary.zlo;

	return boundary;
}

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	FILE *input;
	input = fopen (dataFileName, "r");

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	// DATA_ATOMS *atoms;
	// DATA_BONDS *bonds;
	// DATA_ANGLES *angles;
	// DATA_DIHEDRALS *dihedrals;
	// DATA_IMPROPERS *impropers;
	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);

		if (strstr (lineString, "bond types"))
			sscanf (lineString, "%d \n", &datafile.nBondTypes);

		if (strstr (lineString, "angle types"))
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);

		if (strstr (lineString, "dihedral types"))
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);

		if (strstr (lineString, "improper types"))
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	return datafile;
}

int countNAtoms (FILE *file_dump, int *nAtomEntries)
{
	int nAtoms, currentAtomID, nAtomsFixed;
	char lineString[2000];
	rewind (file_dump);

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, file_dump); }

	sscanf (lineString, "%d\n", &nAtoms);
	(*nAtomEntries) = nAtoms;
	rewind (file_dump);
	nAtomsFixed = nAtoms;

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > nAtoms) {
			nAtomsFixed = currentAtomID; }
	}

	return nAtomsFixed;
}

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtoms, int nAtomEntries, SIMULATION_BOUNDARY *boundary)
{
	char lineString[2000];
	int currentAtomID = 1;

	for (int i = 0; i < 5; ++i) {
		fgets (lineString, 2000, file_dump); }

	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).xlo, &(*boundary).xhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).ylo, &(*boundary).yhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).zlo, &(*boundary).zhi);
	fgets (lineString, 2000, file_dump);

	for (int i = 0; i < nAtomEntries; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);
		sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
		atoms[currentAtomID - 1].isEndGroup = 0;
	}

	return atoms;
}

TRAJECTORY *initializeMolID (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms + 1; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].molType = 0;
	}
	return atoms;
}

TRAJECTORY *assignMolID (TRAJECTORY *atoms, int nAtoms, int *nCTAB, int *nDDAB)
{
	int nAtomsInbetween = 0, IDofPreviousBr = 0;
	(*nCTAB) = 0; (*nDDAB) = 0;
	int nSurfactants = 0;

	for (int i = 0; i < nAtoms; ++i)
	{		
		if (atoms[i].atomType == 5)
		{
			nAtomsInbetween = atoms[i].atomID - IDofPreviousBr;
			IDofPreviousBr = atoms[i].atomID;
			nSurfactants++;

			if (nAtomsInbetween == 63)
			{
				(*nCTAB)++;
				atoms[i].molType = 1;
				for (int j = i; j > (i - 63); --j)
				{
					if (atoms[j].molType == 0) {
						atoms[j].molType = 1; }
				}
			}
			else if (nAtomsInbetween == 84)
			{
				(*nDDAB)++;
				atoms[i].molType = 2;
				for (int j = i; j > (i - 84); --j)
				{
					if (atoms[j].molType == 0) {
						atoms[j].molType = 2; }
				}
			}

			nAtomsInbetween = 0;
		}
	}

	return atoms;
}

TRAJECTORY *countFullVectors (TRAJECTORY *atoms, int nAtoms, DATAFILE_INFO datafile, int minAtomType, int maxAtomType, int *dFullCount)
{
	int nSurfactantAtoms = 0;
	(*dFullCount) = 0;

	/*
		atoms[i].isEndGroup = 1 denotes the tail C atom
		atoms[i].isEndGroup = 2 denotes the head N atom
	*/

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == 2) {
			atoms[i].isEndGroup = 2; }

		if (atoms[i].atomType <= 7) {
			nSurfactantAtoms++; }
	}

	/*
		End tail groups for DDAB are 1 and 27.
		End tail group for CTAB is 20
	*/

	int currentAtom = 0, /*currentMolID = 0,*/ nAtomsInMolecule = 0;

	while (currentAtom < nSurfactantAtoms)
	{
		if (currentAtom == 0)
		{
			if (atoms[0].molType == 1)
			{
				atoms[currentAtom + 19].isEndGroup = 1;
				nAtomsInMolecule = 63;
				currentAtom += nAtomsInMolecule;
				(*dFullCount) += 1;
			}
			if (atoms[0].molType == 2)
			{
				atoms[currentAtom].isEndGroup = 1;
				atoms[currentAtom + 26].isEndGroup = 1;
				nAtomsInMolecule = 84;
				currentAtom += nAtomsInMolecule;
				(*dFullCount) += 1;
			}
		}
		else
		{
			if (atoms[currentAtom].molType == 1)
			{
				atoms[currentAtom + 19].isEndGroup = 1;
				nAtomsInMolecule = 63;
				currentAtom += nAtomsInMolecule;
				(*dFullCount) += 1;
			}
			if (atoms[currentAtom].molType == 2)
			{
				atoms[currentAtom].isEndGroup = 1;
				atoms[currentAtom + 26].isEndGroup = 1;
				nAtomsInMolecule = 84;
				currentAtom += nAtomsInMolecule;
				(*dFullCount) += 2;
			}
		}

		if (atoms[currentAtom].atomType > maxAtomType || atoms[currentAtom].atomType < minAtomType)
		{
			currentAtom++;
		}
	}

	return atoms;
}

// This currently includes H atom. Exclude H atom by some means. Do the same for other vectors that are based on bonds.
// This function will automatically exclude the dihedrals with H atoms when the input lammps trajectory does not contain H atoms
void countD4Vectors (TRAJECTORY *atoms, DATAFILE_INFO datafile, DATA_DIHEDRALS *dataDihedrals, int minAtomType, int maxAtomType, int *d4Count)
{
	(*d4Count) = 0;

	for (int i = 0; i < datafile.nDihedrals; ++i)
	{
		if (atoms[dataDihedrals[i].atom1 - 1].atomType >= minAtomType && atoms[dataDihedrals[i].atom4 - 1].atomType <= maxAtomType) {
			(*d4Count)++; }
	}
}

void countD1Vectors (TRAJECTORY *atoms, DATAFILE_INFO datafile, DATA_BONDS *dataBonds, int minAtomType, int maxAtomType, int *d1Count)
{
	(*d1Count) = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if (atoms[dataBonds[i].atom1 - 1].atomType >= minAtomType && atoms[dataBonds[i].atom2 - 1].atomType <= maxAtomType) {
			(*d1Count)++; }
	}
}

VECTOR *computeFullVectors (VECTOR *dFull, int dFullCount, TRAJECTORY *atoms, int nAtoms)
{
	int currentVector = 0, firstVectorAssigned = 0, secondVectorAssigned = 0;

	for (int i = 0; i < dFullCount; ++i)
	{
		dFull[i].x1 = 0;
		dFull[i].y1 = 0;
		dFull[i].z1 = 0;
		dFull[i].x2 = 0;
		dFull[i].y2 = 0;
		dFull[i].z2 = 0;
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].molType == 1)
		{
			if (firstVectorAssigned == 0 && atoms[i].isEndGroup == 2 && atoms[i].atomType == 2)
			{
				dFull[currentVector].x1 = atoms[i].x;
				dFull[currentVector].y1 = atoms[i].y;
				dFull[currentVector].z1 = atoms[i].z;
				firstVectorAssigned = 1;

				if (currentVector >= dFullCount) {
					goto endThisLoop; }
			}

			if (firstVectorAssigned == 1 && atoms[i].isEndGroup == 1 && atoms[i].atomType == 3)
			{
				dFull[currentVector].x2 = atoms[i].x;
				dFull[currentVector].y2 = atoms[i].y;
				dFull[currentVector].z2 = atoms[i].z;
				currentVector++;
				firstVectorAssigned = 0;
			}
		}
		if (atoms[i].molType == 2)
		{
			if (firstVectorAssigned == 0 && atoms[i].isEndGroup == 2 && atoms[i].atomType == 2)
			{
				dFull[currentVector].x1 = atoms[i].x;
				dFull[currentVector].y1 = atoms[i].y;
				dFull[currentVector].z1 = atoms[i].z;
				firstVectorAssigned = 1;

				if (currentVector >= dFullCount) {
					goto endThisLoop; }
			}
			if (secondVectorAssigned == 0 && atoms[i].isEndGroup == 2 && atoms[i].atomType == 2)
			{
				dFull[currentVector + 1].x1 = atoms[i].x;
				dFull[currentVector + 1].y1 = atoms[i].y;
				dFull[currentVector + 1].z1 = atoms[i].z;
				secondVectorAssigned = 1;

				if (currentVector >= dFullCount) {
					goto endThisLoop; }
			}

			if (firstVectorAssigned == 1 && secondVectorAssigned == 1 && atoms[i].isEndGroup == 1 && atoms[i].atomType == 3)
			{
				dFull[currentVector].x2 = atoms[i].x;
				dFull[currentVector].y2 = atoms[i].y;
				dFull[currentVector].z2 = atoms[i].z;

				dFull[currentVector + 1].x2 = atoms[i].x;
				dFull[currentVector + 1].y2 = atoms[i].y;
				dFull[currentVector + 1].z2 = atoms[i].z;
				currentVector += 2;

				firstVectorAssigned = 0;
				secondVectorAssigned = 0;
			}
		}
	}

	endThisLoop:;

	return dFull;
}

VECTOR *computeD4Vectors (VECTOR *d4, int d4Count, TRAJECTORY *atoms, int nAtoms, DATAFILE_INFO datafile, DATA_DIHEDRALS *dataDihedrals, int minAtomType, int maxAtomType)
{
	int currentVector = 0;

	for (int i = 0; i < datafile.nDihedrals; ++i)
	{
		d4[i].x1 = 0;
		d4[i].y1 = 0;
		d4[i].z1 = 0;
		d4[i].x2 = 0;
		d4[i].y2 = 0;
		d4[i].z2 = 0;
	}

	for (int i = 0; i < datafile.nDihedrals; ++i)
	{
		if (atoms[dataDihedrals[i].atom1 - 1].atomType >= minAtomType && atoms[dataDihedrals[i].atom4 - 1].atomType <= maxAtomType) 
		{
			d4[currentVector].x1 = atoms[dataDihedrals[i].atom1 - 1].x;
			d4[currentVector].y1 = atoms[dataDihedrals[i].atom1 - 1].y;
			d4[currentVector].z1 = atoms[dataDihedrals[i].atom1 - 1].z;

			d4[currentVector].x2 = atoms[dataDihedrals[i].atom4 - 1].x;
			d4[currentVector].y2 = atoms[dataDihedrals[i].atom4 - 1].y;
			d4[currentVector].z2 = atoms[dataDihedrals[i].atom4 - 1].z;

			currentVector++;
		}
	}

	return d4;
}

VECTOR *computeD1Vectors (VECTOR *d1, int d1Count, TRAJECTORY *atoms, int nAtoms, DATAFILE_INFO datafile, DATA_BONDS *dataBonds, int minAtomType, int maxAtomType)
{
	int currentVector = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if (atoms[dataBonds[i].atom1 - 1].atomType >= minAtomType && atoms[dataBonds[i].atom2 - 1].atomType <= maxAtomType) 
		{
			d1[currentVector].x1 = atoms[dataBonds[i].atom1 - 1].x;
			d1[currentVector].y1 = atoms[dataBonds[i].atom1 - 1].y;
			d1[currentVector].z1 = atoms[dataBonds[i].atom1 - 1].z;

			d1[currentVector].x2 = atoms[dataBonds[i].atom2 - 1].x;
			d1[currentVector].y2 = atoms[dataBonds[i].atom2 - 1].y;
			d1[currentVector].z2 = atoms[dataBonds[i].atom2 - 1].z;

			currentVector++;
		}
	}

	return d1;
}

VECTOR *computeVectorCenter (VECTOR *dVector, int dVectorCount, SIMULATION_BOUNDARY boundary)
{

	for (int i = 0; i < dVectorCount; ++i)
	{
		dVector[i].xc = 0; dVector[i].yc = 0; dVector[i].zc = 0;
	}

	for (int i = 0; i < dVectorCount; ++i)
	{
		if (fabs (dVector[i].x1 - dVector[i].x2) > (boundary.xLength / 2))
		{
			if (dVector[i].x1 > dVector[i].x2) {
				dVector[i].x2 += boundary.xLength; }
			else {
				dVector[i].x2 -= boundary.xLength; }
		}

		if (fabs (dVector[i].y1 - dVector[i].y2) > (boundary.yLength / 2))
		{
			if (dVector[i].y1 > dVector[i].y2) {
				dVector[i].y2 += boundary.yLength; }
			else {
				dVector[i].y2 -= boundary.yLength; }
		}

		if (fabs (dVector[i].z1 - dVector[i].z2) > (boundary.zLength / 2))
		{
			if (dVector[i].z1 > dVector[i].z2) {
				dVector[i].z2 += boundary.zLength; }
			else {
				dVector[i].z2 -= boundary.zLength; }
		}

		dVector[i].xc = (dVector[i].x1 + dVector[i].x2) / 2;
		dVector[i].yc = (dVector[i].y1 + dVector[i].y2) / 2;
		dVector[i].zc = (dVector[i].z1 + dVector[i].z2) / 2;

		// TODO: Properly check the minimum image convention
		// TODO: push the bond center within the simulation box bounds, if it is outside
	}

	return dVector;
}

RDF *initializeRDF (RDF *rdfData, int rdf_nBins, float rdf_binSize)
{
	for (int i = 0; i < rdf_nBins; ++i)
	{
		if (i == 0)
		{
			rdfData[i].rlo = 0;
			rdfData[i].rhi = rdfData[i].rlo + rdf_binSize;
		}
		else
		{
			rdfData[i].rlo = rdfData[i - 1].rhi;
			rdfData[i].rhi = rdfData[i].rlo + rdf_binSize;
		}

		rdfData[i].gofr = 0;
	}
	return rdfData;
}

RDF *computeVectorRDF (RDF *rdfData, int rdf_nBins, VECTOR *dVector, int dVectorCount, SIMULATION_BOUNDARY boundary, int maxDistance)
{
	float distance1, distance2, xOffset, yOffset, zOffset;

	// vector1
	for (int i = 0; i < dVectorCount; ++i)
	{
		// vector2
		for (int j = 0; j < dVectorCount; ++j)
		{
			if (i != j)
			{
				distance1 = 0; distance2 = 0; xOffset = 0; yOffset = 0; zOffset = 0;

				while (1)
				{
					if (xOffset == 0 && yOffset == 0 && zOffset == 0)
					{
						distance1 = sqrt (
							pow ((dVector[i].xc - dVector[j].xc), 2) +
							pow ((dVector[i].yc - dVector[j].yc), 2) +
							pow ((dVector[i].zc - dVector[j].zc), 2)
							);

						distance2 = -1;
					}
					else
					{
						distance1 = sqrt (
							pow ((dVector[i].xc - dVector[j].xc - (xOffset * boundary.xLength)), 2) +
							pow ((dVector[i].yc - dVector[j].yc - (yOffset * boundary.yLength)), 2) +
							pow ((dVector[i].zc - dVector[j].zc - (zOffset * boundary.zLength)), 2)
							);

						if (distance1 > maxDistance) {
							distance1 = -1; }

						distance2 = sqrt (
							pow ((dVector[i].xc - dVector[j].xc + (xOffset * boundary.xLength)), 2) +
							pow ((dVector[i].yc - dVector[j].yc + (yOffset * boundary.yLength)), 2) +
							pow ((dVector[i].zc - dVector[j].zc + (zOffset * boundary.zLength)), 2)
							);

						if (distance2 > maxDistance) {
							distance2 = -1; }
					}

					if (distance1 == -1 && distance2 == -1) {
						break; }

					// checking rdf bin
					for (int k = 0; k < rdf_nBins; ++k)
					{
						if (distance1 > rdfData[k].rlo && distance1 <= rdfData[k].rhi) {
							rdfData[k].gofr++; }
						if (distance2 > rdfData[k].rlo && distance2 <= rdfData[k].rhi) {
							rdfData[k].gofr++; }
					}

					xOffset++; yOffset++;
				}
			}
		}
	}

	// for (int i = 0; i < rdf_nBins; ++i)
	// {
	// 	printf("%f %f %f\n", rdfData[i].rlo, rdfData[i].rhi, (rdfData[i].gofr / rdfData[rdf_nBins - 1].gofr * 2));
	// 	usleep (100000);
	// }

	return rdfData;
}

TRAJECTORY *initializeAtoms (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].molType = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].isEndGroup = 0;
	}

	return atoms;
}

float computeCosTheta (float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{
	float dotProduct, magnitude1, magnitude2, cosTheta/*, orderParameter*/;

	dotProduct = (x2 - x1) * (x4 - x3) + (y2 - y1) * (y4 - y3) + (z2 - z1) * (z4 - z3);
	magnitude1 = sqrt (
		(x2 - x1) * (x2 - x1) + 
		(y2 - y1) * (y2 - y1) + 
		(z2 - z1) * (z2 - z1)
		);
	magnitude2 = sqrt (
		(x4 - x3) * (x4 - x3) + 
		(y4 - y3) * (y4 - y3) + 
		(z4 - z3) * (z4 - z3)
		);
	cosTheta = dotProduct / (magnitude1 * magnitude2);

	return cosTheta;
}

STATS computeNormOrderParameter (STATS orderParameter_norm, VECTOR *dFull, int dFullCount, int currentTimeframe, FILE *file_orderParameterNorm)
{
	orderParameter_norm.average = 0;
	orderParameter_norm.standardDeviation = 0;
	float /*dotProduct, magnitude1, magnitude2,*/ cosTheta;
	// float vectorX = 0, vectorY = 0, vectorZ = 1;

	for (int i = 0; i < dFullCount; ++i)
	{
		// dotProduct = (dFull[i].z2 - dFull[i].z1) * (1 - 0);
		// magnitude1 = sqrt (
		// 	(dFull[i].x2 - dFull[i].x1) * (dFull[i].x2 - dFull[i].x1) +
		// 	(dFull[i].y2 - dFull[i].y1) * (dFull[i].y2 - dFull[i].y1) +
		// 	(dFull[i].z2 - dFull[i].z1) * (dFull[i].z2 - dFull[i].z1)
		// 	);
		// magnitude2 = 1;
		// cosTheta = dotProduct / (magnitude1 * magnitude2);
		cosTheta = computeCosTheta (dFull[i].x1, dFull[i].y1, dFull[i].z1, dFull[i].x2, dFull[i].y2, dFull[i].z2, 0, 0, 0, 0, 0, 1);
		orderParameter_norm.average += (cosTheta * cosTheta);
	}

	orderParameter_norm.average /= dFullCount;
	orderParameter_norm.average = 1.5 * orderParameter_norm.average - 0.5;

	fprintf(file_orderParameterNorm, "%d %f\n", currentTimeframe, orderParameter_norm.average);
	fflush (file_orderParameterNorm);

	return orderParameter_norm;
}

float translatePeriodic (float r1, float r2, float simulationBoxLength)
{
	if (fabs (r1 - r2) > (simulationBoxLength / 2))
	{
		if (r1 > r2) {
			r2 += simulationBoxLength; }
		else if (r2 > r1) {
			r2 -= simulationBoxLength; }
	}

	return r2;
}

float computePeriodicDistance (float x1, float y1, float z1, float x2, float y2, float z2, float xLength, float yLength, float zLength)
{
	float distance;

	x2 = translatePeriodic (x1, x2, xLength);
	y2 = translatePeriodic (y1, y2, yLength);
	z2 = translatePeriodic (z1, z2, zLength);

	distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

	return distance;
}

VECTOR *micCorrectionVectors (VECTOR *d, int dCount, SIMULATION_BOUNDARY boundary)
{
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);
	// float previousCoords, translatedCoords;

	for (int i = 0; i < dCount; ++i)
	{
		// previousCoords = dFull[i].x2;
		d[i].x2 = translatePeriodic (d[i].x1, d[i].x2, xLength);
		d[i].y2 = translatePeriodic (d[i].y1, d[i].y2, yLength);
		d[i].z2 = translatePeriodic (d[i].z1, d[i].z2, zLength);
		// translatedCoords = dFull[i].x2;
		// if (previousCoords != translatedCoords)
		// {
		// 	printf("%f -> %f\n", previousCoords, translatedCoords);
		// 	usleep (100000);
		// }
	}

	return d;
}

ORDERPARAMETER_BINS *assignOrderParameterBins (ORDERPARAMETER_BINS *orderParameter_dBins, SIMULATION_BOUNDARY boundary, float orderParameterBinDistance, int *nBins)
{
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);
	float maxLength = sqrt ((xLength * xLength) + (yLength * yLength) + (zLength * zLength));
	(*nBins) = ceil (maxLength / orderParameterBinDistance);

	orderParameter_dBins = (ORDERPARAMETER_BINS *) malloc ((*nBins) * sizeof (ORDERPARAMETER_BINS));

	for (int i = 0; i < (*nBins); ++i)
	{
		if (i == 0)
		{
			orderParameter_dBins[i].rlo = 0;
			orderParameter_dBins[i].rhi = orderParameterBinDistance;
			orderParameter_dBins[i].orderParameter = 0;
			orderParameter_dBins[i].count = 0;
		}
		else
		{
			orderParameter_dBins[i].rlo = orderParameter_dBins[i - 1].rhi;
			orderParameter_dBins[i].rhi = orderParameter_dBins[i].rlo + orderParameterBinDistance;
			orderParameter_dBins[i].orderParameter = 0;
			orderParameter_dBins[i].count = 0;
		}
	}

	return orderParameter_dBins;
}

ORDERPARAMETER_BINS *resetOrderParameterBins (ORDERPARAMETER_BINS *orderParameter, int nBins)
{
	for (int i = 0; i < nBins; ++i)
	{
		orderParameter[i].orderParameter = 0;
		orderParameter[i].count = 0;
	}

	return orderParameter;
}

ORDERPARAMETER_BINS *computeOrderParameterVDistance (ORDERPARAMETER_BINS *orderParameter_dBins, ORDERPARAMETER_BINS *orderParameter_dBins_local, VECTOR *d, int dCount, SIMULATION_BOUNDARY boundary, int nBins_d, int currentTimeframe, FILE *file_dVDistance)
{
	float distance, cosTheta;
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);

	for (int i = 0; i < dCount; ++i)
	{
		if ((i%100) == 0)
		{
			printf("computing order parameter for vector: %d/%d                     \r", i, dCount);
			fflush (stdout);
		}


		for (int j = 0; j < dCount; ++j)
		{
			if (i != j)
			{
				// Make this d[i].x1, d[i].y1, d[i].z1 to d[j].x1, d[j].y1, d[j].z1 and plot with distance
				// x1, y1, z1 belongs to the head groups
				distance = computePeriodicDistance (d[i].xc, d[i].yc, d[i].zc, d[j].xc, d[j].yc, d[j].zc, xLength, yLength, zLength);
				cosTheta = computeCosTheta (d[i].x1, d[i].y1, d[i].z1, d[i].x2, d[i].y2, d[i].z2, d[j].x1, d[j].y1, d[j].z1, d[j].x2, d[j].y2, d[j].z2);

				for (int k = 0; k < nBins_d; ++k)
				{
					if (distance < (float)orderParameter_dBins_local[k].rhi && distance >= (float)orderParameter_dBins_local[k].rlo) {
						orderParameter_dBins_local[k].orderParameter += 1.5 * ((long double)cosTheta * (long double)cosTheta) - 0.5;
						orderParameter_dBins_local[k].count += (long double)1; }
				}
			}
		}
	}

	// Calculating the average cosTheta and then the order parameter

	fprintf(file_dVDistance, "# %d\n", currentTimeframe);

	for (int i = 0; i < nBins_d; ++i)
	{
		orderParameter_dBins[i].orderParameter = (orderParameter_dBins_local[i].orderParameter / orderParameter_dBins_local[i].count);

/*		if (orderParameter_dBins[i].orderParameter != orderParameter_dBins[i].orderParameter) {
			orderParameter_dBins[i].orderParameter = -1; }
		else {
			orderParameter_dBins[i].orderParameter = 1.5 * orderParameter_dBins[i].orderParameter - 0.5; }
*/
		fprintf(file_dVDistance, "%LF %LF %LF\n", orderParameter_dBins[i].rlo, orderParameter_dBins[i].rhi, orderParameter_dBins[i].orderParameter);
	}

	fflush (file_dVDistance);

	return orderParameter_dBins_local;
}

ORDERPARAMETER_BINS *computeOrderParameterVDistance2 (ORDERPARAMETER_BINS *orderParameter_dBins, ORDERPARAMETER_BINS *orderParameter_dBins_local, VECTOR *d, int dCount, SIMULATION_BOUNDARY boundary, int nBins_d, int currentTimeframe, FILE *file_dVDistance)
{
	float distance, cosTheta;
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);

	for (int i = 0; i < dCount; ++i)
	{
		if ((i%100) == 0)
		{
			printf("computing order parameter for vector: %d/%d                     \r", i, dCount);
			fflush (stdout);
		}


		for (int j = 0; j < dCount; ++j)
		{
			if (i != j)
			{
				// Make this d[i].x1, d[i].y1, d[i].z1 to d[j].x1, d[j].y1, d[j].z1 and plot with distance
				// x1, y1, z1 belongs to the head groups
				distance = computePeriodicDistance (d[i].x1, d[i].y1, d[i].z1, d[j].x1, d[j].y1, d[j].z1, xLength, yLength, zLength);
				cosTheta = computeCosTheta (d[i].x1, d[i].y1, d[i].z1, d[i].x2, d[i].y2, d[i].z2, d[j].x1, d[j].y1, d[j].z1, d[j].x2, d[j].y2, d[j].z2);

				for (int k = 0; k < nBins_d; ++k)
				{
					if (distance < (float)orderParameter_dBins_local[k].rhi && distance >= (float)orderParameter_dBins_local[k].rlo) {
						orderParameter_dBins_local[k].orderParameter += ((long double)cosTheta * (long double)cosTheta);
						orderParameter_dBins_local[k].count += (long double)1; }
				}
			}
		}
	}

	// Calculating the average cosTheta and then the order parameter

	fprintf(file_dVDistance, "# %d\n", currentTimeframe);

	for (int i = 0; i < nBins_d; ++i)
	{
		orderParameter_dBins[i].orderParameter = (orderParameter_dBins_local[i].orderParameter / orderParameter_dBins_local[i].count);

		if (orderParameter_dBins[i].orderParameter != orderParameter_dBins[i].orderParameter) {
			orderParameter_dBins[i].orderParameter = -1; }
		else {
			orderParameter_dBins[i].orderParameter = 1.5 * orderParameter_dBins[i].orderParameter - 0.5; }

		fprintf(file_dVDistance, "%LF %LF %LF\n", orderParameter_dBins[i].rlo, orderParameter_dBins[i].rhi, orderParameter_dBins[i].orderParameter);
	}

	fflush (file_dVDistance);

	return orderParameter_dBins_local;
}

VECTOR *computeDDDABVectors (VECTOR *dDDAB, int nDDAB, TRAJECTORY *atoms, int nAtoms)
{
	int currentDDAB = 0;
	bool DDABmolecule;
	DDABmolecule = false;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].molType == 2 && atoms[i].isEndGroup == 1)
		{
			if (DDABmolecule)
			{
				dDDAB[currentDDAB].x2 = atoms[i].x;
				dDDAB[currentDDAB].y2 = atoms[i].y;
				dDDAB[currentDDAB].z2 = atoms[i].z;
				DDABmolecule = false;
				currentDDAB++;
			}
			else
			{
				dDDAB[currentDDAB].x1 = atoms[i].x;
				dDDAB[currentDDAB].y1 = atoms[i].y;
				dDDAB[currentDDAB].z1 = atoms[i].z;
				DDABmolecule = true;
			}
		}
	}

	return dDDAB;
}

int main(int argc, char const *argv[])
{
	FILE *file_dump, *file_data, *file_orderParameterNorm, *file_dFullVDistance, *file_dFullVDistance2, *file_d4VDistance, *file_dDDABVDistance/*, *file_d1VDistance*/;

	char *pipeString;
	pipeString = (char *) malloc (1000 * sizeof (char));

	if (strstr (argv[1], ".xz")) {
		snprintf (pipeString, 1000, "xzcat %s", argv[1]);
		file_dump = popen (pipeString, "r"); }
	else {
		file_dump = fopen (argv[1], "r"); }

	file_data = fopen (argv[2], "r");
	file_orderParameterNorm = fopen ("orderParameter.normal", "w");
	file_dFullVDistance = fopen ("dFullVDistance.orderParameter", "w");
	file_dFullVDistance2 = fopen ("dFullVDistance2.orderParameter", "w");
	file_d4VDistance = fopen ("d4VDistance.orderParameter", "w");
	file_dDDABVDistance = fopen ("dDDABVDistance.orderParameter", "w");
	// file_d1VDistance = fopen ("d1VDistance.orderParameter", "w");

	int file_status, nAtoms, currentTimeframe = 0, nAtomEntries;
	nAtoms = countNAtoms (file_dump, &nAtomEntries);

	if (strstr (argv[1], ".xz")) {
		pclose (file_dump);
		snprintf (pipeString, 1000, "xzcat %s", argv[1]);
		file_dump = popen (pipeString, "r"); }
	else {
		fclose (file_dump);
		file_dump = fopen (argv[1], "r"); }

	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	printf("Number of atoms in the trajectory file: %d\n", nAtoms);

	int dFullCount, d4Count, d1Count;
	int nCTAB, nDDAB;
	VECTOR *dFull, *d4, *d1, *dDDAB;

	RDF *rdf_dFull, *rdf_d4, *rdf_d1;
	rdf_dFull = (RDF *) malloc (100 * sizeof (RDF));
	rdf_d4 = (RDF *) malloc (100 * sizeof (RDF));
	rdf_d1 = (RDF *) malloc (100 * sizeof (RDF));

	rdf_dFull = initializeRDF (rdf_dFull, 100, 3);
	rdf_d4 = initializeRDF (rdf_d4, 100, 3);
	rdf_d1 = initializeRDF (rdf_d1, 100, 3);

	DATA_ATOMS *dataAtoms;
	DATA_BONDS *dataBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	SIMULATION_BOUNDARY boundary;
	boundary = readDumpBoundary (file_dump, boundary);

	if (strstr (argv[1], ".xz")) {
		pclose (file_dump);
		snprintf (pipeString, 1000, "xzcat %s", argv[1]);
		file_dump = popen (pipeString, "r"); }
	else {
		fclose (file_dump);
		file_dump = fopen (argv[1], "r"); }

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);

	atoms = initializeAtoms (atoms, nAtoms);

	STATS orderParameter_norm;

	int nBins_dFull, nBins_d4, nBins_d1, nBins_dDDAB;

	ORDERPARAMETER_BINS *orderParameter_dFullBins, *orderParameter_dFullBins2;
	orderParameter_dFullBins = assignOrderParameterBins (orderParameter_dFullBins, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dFull);
	orderParameter_dFullBins2 = assignOrderParameterBins (orderParameter_dFullBins2, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dFull);

	ORDERPARAMETER_BINS *orderParameter_d4Bins;
	orderParameter_d4Bins = assignOrderParameterBins (orderParameter_d4Bins, boundary, ORDERPARAMETERBINDISTANCE, &nBins_d4);

	ORDERPARAMETER_BINS *orderParameter_d1Bins;
	orderParameter_d1Bins = assignOrderParameterBins (orderParameter_d1Bins, boundary, ORDERPARAMETERBINDISTANCE, &nBins_d1);

	ORDERPARAMETER_BINS *orderParameter_dDDABBins;
	orderParameter_dDDABBins = assignOrderParameterBins (orderParameter_dDDABBins, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dDDAB);

	ORDERPARAMETER_BINS *orderParameter_dFullBins_local, *orderParameter_dFullBins2_local;
	orderParameter_dFullBins_local = assignOrderParameterBins (orderParameter_dFullBins_local, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dFull);
	orderParameter_dFullBins2_local = assignOrderParameterBins (orderParameter_dFullBins2_local, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dFull);

	ORDERPARAMETER_BINS *orderParameter_d4Bins_local;
	orderParameter_d4Bins_local = assignOrderParameterBins (orderParameter_d4Bins_local, boundary, ORDERPARAMETERBINDISTANCE, &nBins_d4);

	ORDERPARAMETER_BINS *orderParameter_d1Bins_local;
	orderParameter_d1Bins_local = assignOrderParameterBins (orderParameter_d1Bins_local, boundary, ORDERPARAMETERBINDISTANCE, &nBins_d1);

	ORDERPARAMETER_BINS *orderParameter_dDDABBins_local;
	orderParameter_dDDABBins_local = assignOrderParameterBins (orderParameter_dDDABBins_local, boundary, ORDERPARAMETERBINDISTANCE, &nBins_dDDAB);

	file_status = 1;

	while (file_status != EOF)
	{
		fprintf(stdout, "computing timestep: %d...                                                            \n", currentTimeframe);
		fflush (stdout);

		atoms = initializeMolID (atoms, nAtoms);
		atoms = readTimestep (file_dump, atoms, nAtoms, nAtomEntries, &boundary);		
		atoms = assignMolID (atoms, nAtoms, &nCTAB, &nDDAB);
		atoms = countFullVectors (atoms, nAtoms, datafile, 1, 5, &dFullCount); // 1 and 5 are the atom type extremes for surfactants

		// Re-assigning the local order parameter bins every timestep
/*		orderParameter_dFullBins_local = resetOrderParameterBins (orderParameter_dFullBins_local, nBins_dFull);
		orderParameter_dFullBins2_local = resetOrderParameterBins (orderParameter_dFullBins2_local, nBins_dFull);

		orderParameter_d4Bins_local = resetOrderParameterBins (orderParameter_d4Bins_local, nBins_d4);
		orderParameter_d1Bins_local = resetOrderParameterBins (orderParameter_d1Bins_local, nBins_d1);

		orderParameter_dDDABBins_local = resetOrderParameterBins (orderParameter_dDDABBins_local, nBins_dDDAB);
*/
		if (currentTimeframe == 0)
		{
			countD4Vectors (atoms, datafile, dataDihedrals, 1, 5, &d4Count);
			countD1Vectors (atoms, datafile, dataBonds, 1, 5, &d1Count);
			
			// printf("Assigning %d mem for dFull array...\n", dFullCount);
			dFull = (VECTOR *) malloc (dFullCount * sizeof (VECTOR));
			d4 = (VECTOR *) malloc (d4Count * sizeof (VECTOR));
			dDDAB = (VECTOR *) malloc (nDDAB * sizeof (VECTOR));
			// d1 = (VECTOR *) malloc (d1Count * sizeof (VECTOR));
		}

		dFull = computeFullVectors (dFull, dFullCount, atoms, nAtoms);
		d4 = computeD4Vectors (d4, d4Count, atoms, nAtoms, datafile, dataDihedrals, 1, 5);
		// d1 = computeD1Vectors (d1, d1Count, atoms, nAtoms, datafile, dataBonds, 1, 5);
		dDDAB = computeDDDABVectors (dDDAB, nDDAB, atoms, nAtoms);

		dFull = computeVectorCenter (dFull, dFullCount, boundary);
		d4 = computeVectorCenter (d4, d4Count, boundary);
		// d1 = computeVectorCenter (d1, d1Count, boundary);
		dDDAB = computeVectorCenter (dDDAB, nDDAB, boundary);

		// MIC correction on vectors
		dFull = micCorrectionVectors (dFull, dFullCount, boundary);
		d4 = micCorrectionVectors (d4, d4Count, boundary);
		// d1 = micCorrectionVectors (d1, d1Count, boundary);
		dDDAB = micCorrectionVectors (dDDAB, nDDAB, boundary);

		// Calculate RDF by replicating the vectors on all three directions.
		// rdf_dFull = computeVectorRDF (rdf_dFull, 100, dFull, dFullCount, boundary, 300);
		// Use the end of first peak as cut-off distance for calculating angle or order parameter.

		// Calculating order parameter between dFull and vector normal to Au surface
		orderParameter_norm = computeNormOrderParameter (orderParameter_norm, dFull, dFullCount, currentTimeframe, file_orderParameterNorm);

		// Calculating order parameter as a function of radial distance
		if (currentTimeframe%2 == 0)
		{
			orderParameter_dFullBins_local = computeOrderParameterVDistance (orderParameter_dFullBins, orderParameter_dFullBins_local, dFull, dFullCount, boundary, nBins_dFull, currentTimeframe, file_dFullVDistance);
			orderParameter_dFullBins2_local = computeOrderParameterVDistance2 (orderParameter_dFullBins2, orderParameter_dFullBins2_local, dFull, dFullCount, boundary, nBins_dFull, currentTimeframe, file_dFullVDistance2);
		}

		if (currentTimeframe%20 == 0)
		{
			orderParameter_d4Bins_local = computeOrderParameterVDistance (orderParameter_d4Bins, orderParameter_d4Bins_local, d4, d4Count, boundary, nBins_d4, currentTimeframe, file_d4VDistance);
		}

		if (currentTimeframe%1 == 0)
		{
			orderParameter_dDDABBins_local = computeOrderParameterVDistance (orderParameter_dDDABBins, orderParameter_dDDABBins_local, dDDAB, nDDAB, boundary, nBins_dDDAB, currentTimeframe, file_dDDABVDistance);
		}
		// orderParameter_d1Bins = computeOrderParameterVDistance (orderParameter_d1Bins, d1, d1Count, boundary, nBins_d1, file_d1VDistance, currentTimeframe);


		// Calculating average order parameter

		file_status = fgetc (file_dump);
		currentTimeframe++;
	}

	fclose (file_data);
	fclose (file_dump);
	fclose (file_orderParameterNorm);
	fclose (file_dFullVDistance);
	fclose (file_d4VDistance);
	fclose (file_dDDABVDistance);

	return 0;
}