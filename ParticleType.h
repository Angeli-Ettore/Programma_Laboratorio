#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include <iostream>

class ParticleType
{
public:
    ParticleType(const char* Name, double Mass, int Charge):
                 fName{Name}, fMass{Mass}, fCharge{Charge} {}

    virtual ~ParticleType()
    {
        delete fName;
    };

    const char* GetName() const;

    double GetMass() const;

    int GetCharge() const;

    virtual double GetWidth() const;

    virtual void Print() const;

private:
    char const* fName;
    double const fMass;
    int const fCharge;
};


#endif
