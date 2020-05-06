#include "EOS.h"

CEOS::CEOS(enEOSTYPE Type,unsigned int EOSID)
:_Type(Type),_EOSID(EOSID)
{
}

CEOS::~CEOS()
{
}

enEOSTYPE CEOS::Type()
{
	return _Type;
}
