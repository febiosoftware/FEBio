

#ifndef KEYGEN_PRODUCT_IDS_H
#define KEYGEN_PRODUCT_IDS_H 1

//
// product_ids.h
//
// Copyright (c) 2010 - University of Utah Software Development Center
//

struct Product_Id_Item {
  const char * name;
  const char * id; // Must be four characters (specifically digits)
};


// Prodcut ID should NEVER be changed.  All new products should just
// be sequentially added to this list.  Make sure to increment the
// NUMBER_OF_PRODUCTS constant.

#define NUMBER_OF_PRODUCTS 1

//                                                         S/W Name                    Product ID        
Product_Id_Item PRODUCT_ID_LIST[NUMBER_OF_PRODUCTS] = {  { "FEBio",                    "0001" } };

#endif // KEYGEN_PRODUCT_IDS_H
