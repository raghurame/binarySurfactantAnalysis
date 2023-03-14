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

int countNAtoms (FILE *file_dump)
{
	int nAtoms;
	char lineString[2000];
	rewind (file_dump);

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}
	sscanf (lineString, "%d\n", &nAtoms);
	rewind (file_dump);

	return nAtoms;
}

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtoms)
{
	char lineString[2000];

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[i].atomID, &atoms[i].atomType, &atoms[i].x, &atoms[i].y, &atoms[i].z, &atoms[i].ix, &atoms[i].iy, &atoms[i].iz);
		atoms[i].isEndGroup = 0;
	}

	return atoms;
}

TRAJECTORY *initializeMolID (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].molType = 0;
	}
	return atoms;
}

TRAJECTORY *assignMolID (TRAJECTORY *atoms, int nAtoms)
{
	int nAtomsInbetween = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		nAtomsInbetween++;

		if (atoms[i].atomType == 5)
		{
			if (nAtomsInbetween == 63)
			{
				atoms[i].molType = 1;
				for (int j = i; j > (i - 63); --j)
				{
					if (atoms[j].molType == 0) {
						atoms[j].molType = 1; }
				}
			}
			else if (nAtomsInbetween == 84)
			{
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

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == 5) {
			atoms[i].isEndGroup = 2; }

		if (atoms[i].atomType <= 7) {
			nSurfactantAtoms++; }
	}

	/*
		End tail groups for DDAB are 1 and 27.
		End tail group for CTAB is 20
	*/

	int currentAtom = 0, currentMolID = 0, nAtomsInMolecule = 0;

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

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].molType == 1)
		{
			if (firstVectorAssigned == 0 && atoms[i].isEndGroup == 1)
			{
				dFull[currentVector].x1 = atoms[i].x;
				dFull[currentVector].y1 = atoms[i].y;
				dFull[currentVector].z1 = atoms[i].z;

				firstVectorAssigned = 1;
				currentVector++;
			}

			if (firstVectorAssigned == 1 && atoms[i].isEndGroup == 2 && atoms[i].atomType == 5)
			{
				dFull[currentVector].x2 = atoms[i].x;
				dFull[currentVector].y2 = atoms[i].y;
				dFull[currentVector].z2 = atoms[i].z;

				firstVectorAssigned = 0;
			}
		}
		if (atoms[i].molType == 2)
		{
			if (firstVectorAssigned == 0 && atoms[i].isEndGroup == 1)
			{
				dFull[currentVector].x1 = atoms[i].x;
				dFull[currentVector].y1 = atoms[i].y;
				dFull[currentVector].z1 = atoms[i].z;

				firstVectorAssigned = 1;
				currentVector++;
			}
			if (secondVectorAssigned == 0 && atoms[i].isEndGroup == 1)
			{
				dFull[currentVector].x1 = atoms[i].x;
				dFull[currentVector].y1 = atoms[i].y;
				dFull[currentVector].z1 = atoms[i].z;

				secondVectorAssigned = 1;
				currentVector++;
			}

			if (firstVectorAssigned == 1 && secondVectorAssigned == 1 && atoms[i].isEndGroup == 2 && atoms[i].atomType == 5)
			{
				dFull[currentVector].x2 = atoms[i].x;
				dFull[currentVector].y2 = atoms[i].y;
				dFull[currentVector].z2 = atoms[i].z;

				firstVectorAssigned = 0;
				secondVectorAssigned = 0;
			}
		}
	}

	return dFull;
}

VECTOR *computeD4Vectors (VECTOR *d4, int d4Count, TRAJECTORY *atoms, int nAtoms, DATAFILE_INFO datafile, DATA_DIHEDRALS *dataDihedrals, int minAtomType, int maxAtomType)
{
	int currentVector = 0;

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

			d1[currentVector].x2 = atoms[dataBonds[i].atom4 - 1].x;
			d1[currentVector].y2 = atoms[dataBonds[i].atom4 - 1].y;
			d1[currentVector].z2 = atoms[dataBonds[i].atom4 - 1].z;

			currentVector++;
		}
	}

	return d1;
}

int main(int argc, char const *argv[])
{
	FILE *file_dump, *file_data;
	file_dump = fopen (argv[1], "r");
	file_data = fopen (argv[2], "r");

	int file_status, nAtoms = countNAtoms (file_dump), currentTimeframe = 0;
	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));

	int dFullCount, d4Count, d1Count;
	VECTOR *dFull, *d4, *d1;

	DATA_ATOMS *dataAtoms;
	DATA_BONDS *dataBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);

	while (file_status != EOF)
	{
		atoms = readTimestep (file_dump, atoms, nAtoms);
		atoms = initializeMolID (atoms, nAtoms);
		atoms = assignMolID (atoms, nAtoms);

		if (currentTimeframe == 0)
		{
			atoms = countFullVectors (atoms, nAtoms, datafile, 1, 5, &dFullCount); // 1 and 5 are the atom type extremes for surfactants
			countD4Vectors (atoms, datafile, dataDihedrals, 1, 5, &d4Count);
			countD1Vectors (atoms, datafile, dataBonds, 1, 5, &d1Count);

			dFull = (VECTOR *) malloc (dFullCount * sizeof (VECTOR));
			d4 = (VECTOR *) malloc (d4Count * sizeof (VECTOR));
			d1 = (VECTOR *) malloc (d1Count * sizeof (VECTOR));
		}

		dFull = computeFullVectors (dFull, dFullCount, atoms, nAtoms);
		d4 = computeD4Vectors (d4, d4Count, atoms, nAtoms, datafile, dataDihedrals, 1, 5);
		d1 = computeD1Vectors (d1, d1Count, atoms, nAtoms, datafile, dataBonds, 1, 5);

		file_status = fgetc (file_dump);
		currentTimeframe++;
	}

	fclose (file_data);
	fclose (file_dump);
	return 0;
}