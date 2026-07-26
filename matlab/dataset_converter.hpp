#pragma once

#include <memory>

#include "mex.h"
#include "spglib.h"

namespace mexutil {

struct DatasetDeleter {
    void operator()(SpglibDataset* dataset) const noexcept;
};

struct MagneticDatasetDeleter {
    void operator()(SpglibMagneticDataset* dataset) const noexcept;
};

using DatasetPtr = std::unique_ptr<SpglibDataset, DatasetDeleter>;
using MagneticDatasetPtr =
    std::unique_ptr<SpglibMagneticDataset, MagneticDatasetDeleter>;

mxArray* makeDatasetStruct(SpglibDataset const& dataset);
mxArray* makeMagneticDatasetStruct(SpglibMagneticDataset const& dataset);

}  // namespace mexutil
