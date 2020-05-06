#include "Section.h"

CSection::CSection(enSECTIONTYPE Type,unsigned int SECID)
:_Type(Type),_SECID(SECID)
{
}

CSection::~CSection()
{
}

enSECTIONTYPE CSection::Type()
{
	return _Type;
}
