#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

#include <cstdio>

#define allocation_tableau(nom, type, nombre) \
if( (nom = (type*) calloc (nombre, sizeof(type) ) ) == NULL ) \
{\
 printf("\n Allocation dynamique impossible pour un pointeur-tableau \n");\
 exit(1);\
}

void ignorer_commentaires(_IO_FILE* file);
void ecrire_image_ppm(char* filename, unsigned char* data, int width, int height);
void lire_nb_lignes_colonnes_image_ppm(char* filename, int* width, int* height);
void lire_image_ppm(char* filename, unsigned char* data, int size);
void planR(unsigned char* src, unsigned char* dest, int size);
void planV(unsigned char* src, unsigned char* dest, int size);
void planB(unsigned char* src, unsigned char* dest, int size);
void ecrire_image_pgm(char* filename, unsigned char* data, int width, int height);
void lire_nb_lignes_colonnes_image_pgm(char* filename, int* width, int* height);
void lire_image_pgm(char* filename, unsigned char* data, int size);
typedef unsigned char OCTET;
#endif // IMAGE_UTILS_H
