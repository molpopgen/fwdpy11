#pragma once

#include <exception>
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

struct ExtremeMigrationUntilOneGenerationAgo
{
    std::string yaml;
    ExtremeMigrationUntilOneGenerationAgo();
};

struct BadEpochRounding02
{
    std::string yaml;
    BadEpochRounding02();
};

struct NonIntegerStartSize
{
    std::string yaml;
    NonIntegerStartSize();
};

struct NonIntegerEndSize
{
    std::string yaml;
    NonIntegerEndSize();
};

struct LinearSizeChange
{
    std::string yaml;
    LinearSizeChange();
};

struct DemeSizeIsOne
{
    std::string yaml;
    DemeSizeIsOne();
};
