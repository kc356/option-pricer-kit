#include "pricing/MCSettings.h"

namespace OptionPricer
{

    MCSettings::MCSettings(int n_paths, int n_steps, unsigned int seed, bool antithetic)
        : n_paths_(n_paths), n_steps_(n_steps), seed_(seed), antithetic_(antithetic)
    {
    }

} // namespace OptionPricer
