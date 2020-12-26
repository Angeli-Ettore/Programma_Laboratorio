#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H

#include "ParticleType.h"

class ResonanceType : public ParticleType
{
public:
    ResonanceType(const char* Name, double Mass, int Charge, double Width) :
                  ParticleType(Name, Mass, Charge), fWidth{Width} {}

    double GetWidth() const override;

    void Print() const override;

private:
    double const fWidth;
};

#endif
