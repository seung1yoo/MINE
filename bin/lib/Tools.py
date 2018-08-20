#!/usr/bin/python

import os
import json


class App :
    def __init__(self,eco,name) :
        self.tools_home=eco["tools"]["meta"]["path"]["tools"]
        self.env=eco["tools"][name]
        self.home=self.env["HOME"].replace("[TOOLS_HOME]",self.tools_home)
        self.exe=self.env["EXEC"].replace("[HOME]",self.home)
        self.n_cpu=1
        self.que=None
        self.temp_path=None
        self.interpreter=None
        #
        self.n_cpu_s={}
        self.name=self.env["NAME"]
        #
        self.exec_s={}
        if "EXEC_S" in self.env :
            for key in self.env["EXEC_S"] :
                self.exec_s[key]=self.env["EXEC_S"][key].replace("[HOME]",self.home)

        #
        self.version=None
        if "VERSION" in self.env :
            self.version=self.env["VERSION"]
        #
        if "QUE" in self.env :
            self.que=self.env["QUE"]
        #
        if "N_CPU" in self.env :
            self.n_cpu=int(self.env["N_CPU"])
        if "N_CPU_S" in self.env :
            self.n_cpu_s=self.env["N_CPU_S"]
        #
        if "RAM" in self.env :
            self.ram=self.env["RAM"]
        #
        self.profile=None
        if "PROFILE" in self.env :
            self.profile=self.env["PROFILE"]
        #
        self.citation=None
        if "CITATION" in self.env :
            self.citation=self.env["CITATION"]
        #
        if "TEMP_PATH" in self.env :
            self.temp_path=self.env["TEMP_PATH"]
            self.temp_path=self.temp_path.replace("[PID]","%d"%(os.getpid()))
            self.temp_path=self.temp_path.replace("[NAME]",self.name)
        #
        self.interpreter=None
        if "INTERPRETER" in self.env :
            self.interpreter=self.env["INTERPRETER"]

    def to_path(self,path) :
        return path.replace("[HOME]",self.home)

    def get_n_cpu(self,profile=None) :
        if profile==None :
            return self.n_cpu
        return self.n_cpu_s[profile]

    def get_que(self):
        return self.que

    def get_temp_path(self,tag=None,mid=None) :
        temp_path=self.temp_path
        if mid!=None :
             temp_path=temp_path.replace("[MID]",str(mid))
        if tag!=None :
             temp_path="%s/%s"%(self.temp_path,tag)
        return temp_path

    def java_prefix(self,tag=None):
        pre=[self.interpreter,"-Xmx%s"%self.ram]
        if self.temp_path!=None :
            temp_path=self.temp_path
            if tag!=None :
                 temp_path="%s/%s"%(self.temp_path,tag)
            pre.append("-Djava.io.tmpdir=%s"%(temp_path))
        pre.extend(["-jar",self.exe])
        return pre

    def get_profile(self,key) :
        return self.profile[key]


