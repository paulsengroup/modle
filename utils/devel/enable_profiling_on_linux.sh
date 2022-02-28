#!/bin/sh

sudo sh -c 'echo 1 > /proc/sys/kernel/perf_event_paranoid && echo 0 > /proc/sys/kernel/kptr_restrict'
