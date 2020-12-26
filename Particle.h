#ifndef PARTICLE_H
#define PARTICLE_H

#include "ResonanceType.h"

#include <cmath>
#include <cstdlib>

class Particle
{
public:
    Particle(const char* Name, double Px, double Py, double Pz) : fPx{Px}, fPy{Py}, fPz{Pz}
    {
        fIParticle = FindParticle(Name);
        if(fIParticle == fNParticleType)
            std::cout << "Particle " << Name << " not on the list, please add it" << '\n';
    }
    Particle() = default;

    int GetIndex() const;
    double GetImpulseX() const;
    double GetImpulseY() const;
    double GetImpulseZ() const;
    double GetParticleMass() const;
    double GetParticleImpulse() const;
    double GetParticleEnergy() const;
    double GetInvMass(Particle const&p) const;
    double GetInvMass(Particle *p) const;
    const char* GetParticleName() const;
    int GetParticleCharge() const;


    void SetParticleID(int index);
    void SetParticleName(char* name);
    void SetImpulse(double Px, double Py, double Pz);

    void PrintParticle() const;
    int Decay2body(Particle &dau1, Particle &dau2) const;
    
    static void PrintArray();
    static void AddParticle(const char* Name, double Mass, int Charge, double Width = 0.);

private:
    static const int fMaxNumParticleType = 10;
    static ParticleType *fParticleType[fMaxNumParticleType];
    static int fNParticleType;
    int fIParticle = 0;
    double fPx = 0;
    double fPy = 0;
    double fPz = 0;
    static int FindParticle(const char* Name);
    void Boost(double Bx, double By, double Bz);
};

#endif
