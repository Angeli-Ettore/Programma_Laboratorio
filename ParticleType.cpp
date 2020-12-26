#include "ParticleType.h"


const char* ParticleType::GetName() const
{
    return fName;
}
double ParticleType::GetMass() const
{
    return fMass;
}
int ParticleType::GetCharge() const
{
    return fMass;
}
double ParticleType::GetWidth() const
{
    return 0.;
}
void ParticleType::Print() const
{
    std::cout << "Name: " << fName << '\n';
    std::cout << "Mass: " << fMass << '\n';
    std::cout << "Charge: " << fCharge << '\n';
}


