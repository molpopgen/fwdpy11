#pragma once

#include <string>

struct SingleDemeModel
{
    std::string yaml;
    SingleDemeModel();
};

struct SingleDemeModelOneSizeChange
{
    std::string yaml;
    SingleDemeModelOneSizeChange();
};

struct TwoDemePerpetualIslandModel
{
    std::string yaml;
    TwoDemePerpetualIslandModel();
};

struct TwoDemePerpetualIslandModelWithSizeChangeAndExtinction
{
    std::string yaml;
    TwoDemePerpetualIslandModelWithSizeChangeAndExtinction();
};

struct TwoDemesUnequalMerge
{
    std::string yaml;
    TwoDemesUnequalMerge();
};

struct VeryRecentPulseTwoGenerationsAgo
{
    std::string yaml;
    VeryRecentPulseTwoGenerationsAgo();
};

struct ExtremeMigrationUntilOneGenerationAgo {
    std::string yaml;
    ExtremeMigrationUntilOneGenerationAgo();
};
