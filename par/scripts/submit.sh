#!/bin/sh

ssh -l koehn -f -T -L 20000:tea.geophysik.tu-freiberg.de:22 koehn@sshproxy.tu-freiberg.de sleep 300

