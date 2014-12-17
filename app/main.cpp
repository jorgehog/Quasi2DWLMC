#include <iostream>

#include "../WLMC/include/WLMC.h"

#include "../../../Dropbox/kMC_dag/find_energy_function/statedensity.h"

using namespace WLMC;
using namespace arma;
using namespace std;

class Quasi2DSystem : public System
{
    const uint m_L;
    const uint m_sMax;

    ivec m_heights;
    ivec m_trialHeights;
    imat m_presetHeights;

    // System interface
public:
    Quasi2DSystem(const uint L, const uint sMax) :
        //nParticles, NX, NY, NZ, samplesPerCheck, flatnesscriterion, adaptive window stuff.., rng
        System(L, L, 1, 1, 20000, 0.85, 50, 100, 0.1, -1000000000, ".", [] () {return as_scalar(randu(1,1));}, false, false),
        m_L(L),
        m_sMax(sMax),
        m_heights(L, fill::zeros),
        m_trialHeights(L, fill::zeros),
        m_presetHeights(L, 0)
    {

    }

    uint nextSite(const uint particleIndex) const
    {
        return (particleIndex + m_heights.size() + 1)%m_heights.size();
    }

    uint prevSite(const uint particleIndex) const
    {
        return (particleIndex + m_heights.size() - 1)%m_heights.size();
    }

    bool isOccupiedLoction(const uint x, const uint y, const uint z) const
    {
        (void) x;
        (void) y;
        (void) z;

        return false;
    }

    double heightDifference(const ivec &heights, const uint particleIndex) const
    {
        return abs(heights(particleIndex) - heights(nextSite(particleIndex)));
    }

    double getValue(const uint particleIndex) const
    {
        return 0.5*heightDifference(m_heights, particleIndex);
    }

    void onSuggestedTrialMove(const uint particleIndex, const uint xd, const uint yd, const uint zd)
    {
        (void) particleIndex;
        (void) xd;
        (void) yd;
        (void) zd;

        m_trialHeights = m_heights;

        uint nTrials = 10;
        uint particle, x, y, z;
        for (uint i = 0; i < nTrials; ++i)
        {
            getRandomParticleAndDestination(particle, x, y, z);

            bool up = URNG() < 0.5;

            if (up)
            {
                m_trialHeights(particle)++;
            }
            else
            {
                m_trialHeights(particle)--;
            }

        }
    }

    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd)
    {
        (void) particleIndex;
        (void) xd;
        (void) yd;
        (void) zd;

        m_heights = m_trialHeights;
    }

    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const
    {
        x = particleIndex;
        y = 1;
        z = 1;
    }

    double setMaximumConfiguration()
    {
        uint L = m_heights.size();

        m_heights(0) = m_heights(L - 1) = 0;

        int dh = m_sMax/(L-2);
        uint spare = m_sMax - (L-2)*dh;

        for (uint i = 1; i < L - 1; ++i)
        {
            m_heights(i) = (-1 + 2*(i%2))*dh;
        }

        for (uint i = 0; i < spare; ++i)
        {
            if (i%2 == 0)
            {
                m_heights(i)--;
            }
            else
            {
                m_heights(i)++;
            }
        }

        return m_sMax;
    }

    double setMinimumConfiguration()
    {
        m_heights.zeros();

        return 0;
    }

    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd)
    {
        (void) particleIndex;
        (void) xd;
        (void) yd;
        (void) zd;

        double prev = getTotalValue();

        double trialValue = 0;

        for (uint i = 0; i < nParticles(); ++i)
        {
            trialValue += 0.5*heightDifference(m_trialHeights, i);
        }

        double diff = trialValue - prev;


        return diff;

    }

    void savePositionData(const uint framenumber) const
    {
        (void) framenumber;
        m_heights.save("heighmap.arma");
    }

    void onPresetSave(const uint n)
    {
        m_presetHeights.col(n) = m_heights;
    }

    void onPresetLoad(const uint n)
    {
        m_heights = m_presetHeights.col(n);
    }

    void init(const uint nbins)
    {
        m_presetHeights.set_size(m_heights.size(), numberOfPresets(nbins));
    }

    void shiftDOS(vec &logDOS)
    {
        logDOS += log(m_L*(m_L-1)) - logDOS(1);
    }

};

int main()
{

    auto seed = time(NULL);
//        auto seed = 1418837610;
    cout << seed << endl;
    srand(seed);

    uint L = 10;
    //    uint sMax = 17*m_heights.size();
    uint sMax = 100;

    EspensAwesomeFormula<uint> f(L, sMax);
    f.initialize();
    vec aws(sMax+1);

    for (uint s = 0; s <= sMax; ++s)
    {
        aws(s) = f.logMT<double>(s);
    }

    aws.save("awesome.arma");

    Quasi2DSystem system(L, sMax);

    uint nbins = 100;
    bool adaptiveWindows = false;

    double fInit = 1;
    double fEnd = 0.00001;

    system.init(nbins);
    system.execute(nbins, adaptiveWindows, fInit, fEnd);

    return 0;
}

