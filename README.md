# Monet
A Matlab library for studying real space embryo development

The community is producing an ever-greater number of methods for analyzing single cell data. Such data has the format that an entry represents a snapshot of a cell’s state at given time. The purpose of single cell analysis is to look at a collection such snapshots and, by looking at the relative differences between them, learn things about the underlying developmental process. These methods have been very successful, but suffer from the limitation that they aren’t compatible with spatial data – the choreography of where cells go over time – because it is of a fundamentally different format.

We solved this problem by creating a method for converting imaging data depicting the real space development of an embryo into a feature representation analogous in form to single cell data, making it possible to analyze morphogenesis using single cell techniques. The essence of the method is to define a feature representation of a cell’s spatiotemporal context at some point in time, then amass representations of all the cell in an embryo into a dataset and analyze it. That feature representation can be computed using a movie of an embryo’s development over some time period in which every cell is represented as a point, or it can include tracking and lineage information, too. In either case, our method enables scalable, automated phenotyping of embryo morphology and development.  

This document describes our process for designing and benchmarking our method, using C. elegans as a model. 
