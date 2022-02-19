# Yurt-Sketch, a framework for weight scalable graph stream summarization

## Introduction

Graph stream, which refers to a dynamic graph model constructed
by an unbounded edge sequence arriving quickly over time, playing a important role in many scenarios like cyber security, social
network and cloud log analysis. As the graph stream has the characteristics of vast volume, fast update speed, traditional graph storage
models like the adjacency matrix and the adjacency list are no
longer efficiency. Prior art of graph stream summarization like CM,
TCM, GSS etc. use a sketch to compress the graph stream and offer query results with bounded errors. However, all these sketch
allocate a non-scalable size of counter to summarize the weight of
the edge which will face the problem of counter overflow as the
the volume of the graph stream is not always predictable and may
fluctuate over time. To address this issue, we propose Yurt sketch, a weight
scalable sketch with memory efficiency for graph stream summarization especially for skewed datasets. Yurt sketch implements a weight scalable graph stream summarization and reduces
the ARE (average relative error) by a ratio of 70% and 33% for different queries on average.

## Detailed design

The files DESCRIPTION.docx shows the detailed design of Yurt

## Source code

We have implemented Yurt on TCM and GSS in C++. We complete the code on Linux 5.4.0-99-generic and compile successfully using gcc 7.5.0.
