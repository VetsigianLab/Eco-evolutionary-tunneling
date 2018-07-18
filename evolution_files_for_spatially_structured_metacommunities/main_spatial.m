[simulation]=generateInput_Spatial;

simulation = evolve_spatial('result', simulation);

[time]=FindStableCommunity_spatial( simulation );
