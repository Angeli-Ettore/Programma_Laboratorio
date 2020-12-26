#include "Particle.h"

#include <cstdlib>

int Particle::fNParticleType = 0;
ParticleType* Particle::fParticleType[fMaxNumParticleType];

int Particle::GetIndex() const
{
    return fIParticle;
}
double Particle::GetImpulseX() const
{
    return fPx;
}
double Particle::GetImpulseY() const
{
    return fPy;
}
double Particle::GetImpulseZ() const
{
    return fPz;
}
double Particle::GetParticleMass() const
{
    return fParticleType[fIParticle] -> GetMass();
}
double Particle::GetParticleImpulse() const
{
    return sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
}
double Particle::GetParticleEnergy() const
{
    return sqrt(pow(GetParticleMass(), 2) + pow(GetParticleImpulse(), 2));
}
double Particle::GetInvMass(Particle const&p) const
{
    double E = pow(p.GetParticleEnergy() + GetParticleEnergy(), 2);
    double P = pow(GetImpulseX() + p.GetImpulseX(), 2) + pow(GetImpulseY() + p.GetImpulseY(), 2) + pow(GetImpulseZ() + p.GetImpulseZ(), 2);
    return sqrt(E - P);
}
double Particle::GetInvMass(Particle *p) const
{
    double E = pow(p -> GetParticleEnergy() + GetParticleEnergy(), 2);
    double P = pow(GetImpulseX() + p -> GetImpulseX(), 2) + pow(GetImpulseY() + p -> GetImpulseY(), 2) + pow(GetImpulseZ() + p -> GetImpulseZ(), 2);
    return sqrt(E - P);
}
const char* Particle::GetParticleName() const
{
    return fParticleType[fIParticle] -> GetName();
}
int Particle::GetParticleCharge() const
{
    return fParticleType[fIParticle] -> GetCharge();
}



void Particle::SetParticleID(int Index)
{
    if(Index >= fNParticleType)
        std::cout << "Out of range index" << '\n';
    fIParticle = Index;
}
void Particle::SetParticleName(char* Name)
{
    int Index = FindParticle(Name);
    if(Index == fNParticleType)
        std::cout << "Particle name is not in the pre-defined array" << '\n';
    else 
        fIParticle = Index;
}
void Particle::SetImpulse(double Px, double Py, double Pz)
{
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}


void Particle::PrintParticle() const
{
    std::cout << "Particle index: " << GetIndex() << '\n';
    std::cout << "Particle name: " << fParticleType[GetIndex()]->GetName() << '\n';
    std::cout << "Particle impulse: " << GetParticleImpulse() << '\n';
    std::cout << "X component: " << GetImpulseX() << '\n';
    std::cout << "Y component: " << GetImpulseY() << '\n';
    std::cout << "Z component: " << GetImpulseZ() << '\n';
}
int Particle::Decay2body(Particle &dau1,Particle &dau2) const 
{
    if(GetParticleMass() == 0.0)
    {
        std::cout << "Decayment cannot be preformed if mass is zero" << '\n';
        return 1;
    }
    double massMot = GetParticleMass();
    double massDau1 = dau1.GetParticleMass();
    double massDau2 = dau2.GetParticleMass();

    if(fIParticle > -1)
    {
        float x1, x2, w, y1, y2;
        
        double invnum = 1.0/RAND_MAX;
        do 
        {
            x1 = 2.0 * rand()*invnum - 1.0;
            x2 = 2.0 * rand()*invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );
        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        massMot += fParticleType[fIParticle]->GetWidth() * y1;
    }

    if(massMot < massDau1 + massDau2)
    {
        std::cout << "Decayment cannot be preformed because mass is too low in this channel" << '\n';
        std::cout << "Child1 mass: " << massDau1 << '\n';
        std::cout << "Child2 mass: " << massDau2 << '\n';
        return 2;
    }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.SetImpulse(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetImpulse(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}

void Particle::PrintArray()
{
    std::cout << "This array contains " << fNParticleType << " particle." << '\n';
    for(int i = 0; i != fNParticleType; ++i)
    {
        std::cout << "Particle number " << (i + 1) << ":" << '\n';
        fParticleType[i]->Print();
        std::cout << '\n';
    }
}
void Particle::AddParticle(const char* Name, double Mass, int Charge, double Width)
{
    if(Width == 0.)
    {
        if(FindParticle(Name) != fNParticleType)
        {
            std::cout << "Error, "<< Name << " already in the array" << '\n';
        }
        else if(fNParticleType == fMaxNumParticleType)
            std::cout << "There is no space left in the array" << '\n';
        else
        {
            fParticleType[fNParticleType] = new ParticleType{Name, Mass, Charge};
            ++fNParticleType;
            std::cout << "Particle " << Name << " added" << '\n';
            std::cout << "Current number of particle: " << fNParticleType << '\n';
        }
    }
    else
    {
        if(FindParticle(Name) != fNParticleType)
        {
            std::cout << "Error, "<< Name << " already in the array" << '\n';
        }
        else if(fNParticleType == fMaxNumParticleType)
            std::cout << "There is no space left in the array" << '\n';
        else
        {
            fParticleType[fNParticleType] = new ResonanceType{Name, Mass, Charge, Width};
            ++fNParticleType;
            std::cout << "Particle " << Name << " added" << '\n';
            std::cout << "Current number of particle: " << fNParticleType << '\n';
        }
    }
}


int Particle::FindParticle(const char* Name)
{
    for(int i = 0; i != fNParticleType; ++i)
    {
        if(fParticleType[i]->GetName() == Name)
            return i;
    }
    return fNParticleType;
}
void Particle::Boost(double Bx, double By, double Bz)
{
    double energy = GetParticleEnergy();
    double b2 = Bx * Bx + By * By + Bz * Bz;
    double gamma = (1.0 / sqrt(1.0 - b2));
    double bp = Bx * fPx + By * fPy + Bz * fPz;
    double gamma2 = (b2 > 0) ? (gamma - 1.0) / b2
                             : 0.0;
    fPx += (gamma2 * bp * Bx + gamma * Bx * energy);
    fPy += (gamma2 * bp * By + gamma * By * energy);
    fPz += (gamma2 * bp * Bz + gamma * Bz * energy);
}









