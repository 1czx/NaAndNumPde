#pragma once 
#include "TimeIntegrator.h"
#include <map>

class TimeIntegratorFactory{

    public:
    using CreateTimeIntegratorCallBack = TimeIntegrator* (*)();
    
    private:
    using CallbackMap = std::map<std::string, CreateTimeIntegratorCallBack>;
    
    public:
    
    void registerTimeIntegrator(const std::string &ID, CreateTimeIntegratorCallBack createFn){
        callbacks_[ID] = createFn;
    }
    TimeIntegrator* createTimeIntegrator(const std::string &ID){
        if(!callbacks_.count(ID)){
            std::cerr << "TimeIntegratorFactory:: No such TimeIntegrator called '" << ID << "'." << endl;
            return nullptr;
        }
        return callbacks_[ID]();
    }
    
    private:

    CallbackMap callbacks_;
    
    TimeIntegratorFactory() = default;
    TimeIntegratorFactory(const TimeIntegratorFactory&) = default;
    TimeIntegratorFactory& operator = (const TimeIntegratorFactory&) = default;
    ~TimeIntegratorFactory() = default;
    
    public:
    static TimeIntegratorFactory& getInstance(){
        static TimeIntegratorFactory singleton;
        return singleton;
    }
};

static void registerABMs(void)__attribute__((constructor));

void registerABMs(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("ABMs", [](){ return (TimeIntegrator*) new ABMs(); });
}

static void registerAMMs(void)__attribute__((constructor));

void registerAMMs(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("AMMs", [](){ return (TimeIntegrator*) new AMMs(); });
}

static void registerBDFs(void)__attribute__((constructor));

void registerBDFs(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("BDFs", [](){ return (TimeIntegrator*) new BDFs(); });
}

static void registerClassicalRK(void)__attribute__((constructor));

void registerClassicalRK(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("classicalRK", [](){ return (TimeIntegrator*) new classicalRK(); });
}

static void registerESDIRK(void)__attribute__((constructor));

void registerESDIRK(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("ESDIRK", [](){ return (TimeIntegrator*) new ESDIRK(); });
}

static void registerGaussLegendreRKMs(void)__attribute__((constructor));

void registerGaussLegendreRKMs(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("GaussLegendreRKMs", [](){ return (TimeIntegrator*) new GaussLegendreRKMs(); });
}

static void registerFehlberg45(void)__attribute__((constructor));

void registerFehlberg45(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("Fehlberg45", [](){ return (TimeIntegrator*) new Fehlberg45(); });
}

static void registerDormandPrince54(void)__attribute__((constructor));

void registerDormandPrince54(){
    auto& factory = TimeIntegratorFactory::getInstance();
    factory.registerTimeIntegrator("DormandPrince54", [](){ return (TimeIntegrator*) new DormandPrince54(); });
}