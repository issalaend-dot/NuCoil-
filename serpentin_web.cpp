#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <string>
#include <map>
#include <csignal>
#include <cstring>
#include <cstdlib>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>

constexpr double PI = 3.14159265358979323846;

struct ParametrosProceso {
    int modoEntalpia = 1;            // 1: Estequiometrico (combustiones), 2: Directo
    double dhCombSustrato = 2805.0;  // kJ/mol sustrato (glucosa)
    double dhCombBiomasa = 560.0;    // kJ/C-mol biomasa
    double rendimientoYxs = 3.0;     // C-mol biomasa / mol sustrato
    double dhDirecto = 460.0;        // kJ/mol referencia
    double tasaMolar = 15.0;         // mol/h
    double factorCalor = 1.0;        // Fraccion metabolica liberada

    double diametroTanque = 1.20;    // m
    double volumenLiquido = 1.10;    // m^3
    double diametroImpulsor = 0.40;  // m
    double velocidadAgitacion = 2.5; // rps (150 rpm)
    double numeroPotencia = 5.0;     // Np (Rushton)
    double flujoAire = 0.010;        // m^3/s
    double rhoCaldo = 1020.0;        // kg/m^3
    double muCaldo = 0.0020;         // Pa*s
    double kCaldo = 0.60;            // W/(m*K)
    double cpCaldo = 4100.0;         // J/(kg*K)

    double tCaldo = 30.0;            // °C
    double tRefIn = 10.0;            // °C
    double tRefOut = 18.0;           // °C
    double flujoRefDisp = 0.25;      // kg/s
    double cpRef = 4185.0;           // J/(kg*K)
    double rhoRef = 999.0;           // kg/m^3
    double muRef = 0.00117;          // Pa*s
    double kRef = 0.589;             // W/(m*K)

    double dExtTubo = 0.025;         // m
    double espesorPared = 0.002;     // m
    double kPared = 16.0;            // W/(m*K)
    double factorEnsuciamiento = 0.0002; // m^2*K/W
    double holguraPared = 0.10;      // m
    double pasoVueltas = 0.075;      // m
};

struct ResultadosDiseno {
    bool valido = true;
    std::string mensajeError;

    double dhReaccion = 0.0;
    double qReaccion = 0.0;
    double reAgitacion = 0.0;
    double potenciaNoAir = 0.0;
    double potenciaAireada = 0.0;
    double qTotal = 0.0;

    double dtlm = 0.0;
    double ho = 0.0;
    double flujoRefReq = 0.0;
    double velRef = 0.0;
    double reRef = 0.0;
    double hi = 0.0;
    double uGlobal = 0.0;

    double areaReq = 0.0;
    double longitudTubo = 0.0;
    double diametroHelice = 0.0;
    double numVueltas = 0.0;
    double alturaSerpentin = 0.0;
    double alturaLiquido = 0.0;
    double alturaTanqueEst = 0.0;
    double caidaPresion = 0.0;
    bool cabeEnLiquido = true;
    bool libreImpulsor = true;
    bool caudalSuficiente = true;
};

ResultadosDiseno calcularDiseno(const ParametrosProceso& p) {
    ResultadosDiseno r;

    if (p.diametroTanque <= 0.05 || p.volumenLiquido <= 0.001 || p.diametroImpulsor <= 0.01 ||
        p.velocidadAgitacion <= 0.0 || p.dExtTubo <= 0.002 || p.pasoVueltas <= 0.002 ||
        p.rhoCaldo <= 0.0 || p.muCaldo <= 0.0 || p.kCaldo <= 0.0 || p.cpCaldo <= 0.0 ||
        p.rhoRef <= 0.0 || p.muRef <= 0.0 || p.kRef <= 0.0 || p.cpRef <= 0.0 || p.kPared <= 0.0) {
        r.valido = false;
        r.mensajeError = "Todas las dimensiones y propiedades fisicas deben ser mayores a cero.";
        return r;
    }

    if (p.tCaldo <= p.tRefIn || p.tCaldo <= p.tRefOut || p.tRefOut <= p.tRefIn) {
        r.valido = false;
        r.mensajeError = "Error termodinamico: Se requiere T_caldo > T_salida > T_entrada.";
        return r;
    }

    double dIntTubo = p.dExtTubo - 2.0 * p.espesorPared;
    if (dIntTubo <= 0.001) {
        r.valido = false;
        r.mensajeError = "Error geometrico: El espesor de pared supera el radio interior del tubo.";
        return r;
    }

    // 1. Balance de energia y entalpia de reaccion (Doran Cap. 5 y 8)
    if (p.modoEntalpia == 1) {
        r.dhReaccion = p.dhCombSustrato - (p.rendimientoYxs * p.dhCombBiomasa);
    } else {
        r.dhReaccion = p.dhDirecto;
    }
    if (r.dhReaccion < 0.0) r.dhReaccion = 0.0;

    r.qReaccion = (r.dhReaccion * 1000.0 * p.tasaMolar * p.factorCalor) / 3600.0;
    r.reAgitacion = (p.rhoCaldo * p.velocidadAgitacion * std::pow(p.diametroImpulsor, 2)) / p.muCaldo;
    r.potenciaNoAir = p.numeroPotencia * p.rhoCaldo * std::pow(p.velocidadAgitacion, 3) * std::pow(p.diametroImpulsor, 5);

    if (p.flujoAire > 0.0) {
        double term = (std::pow(r.potenciaNoAir, 2) * p.velocidadAgitacion * std::pow(p.diametroImpulsor, 3))
                      / std::pow(p.flujoAire, 0.56);
        r.potenciaAireada = 0.706 * std::pow(term, 0.45);
        if (r.potenciaAireada > r.potenciaNoAir) r.potenciaAireada = r.potenciaNoAir;
    } else {
        r.potenciaAireada = r.potenciaNoAir;
    }

    r.qTotal = r.qReaccion + r.potenciaAireada;
    if (r.qTotal <= 1.0) {
        r.valido = false;
        r.mensajeError = "La carga termica total debe ser positiva para dimensionar el enfriamiento.";
        return r;
    }

    // 2. Fuerza impulsora DTLM
    double dt1 = p.tCaldo - p.tRefIn;
    double dt2 = p.tCaldo - p.tRefOut;
    r.dtlm = (std::abs(dt1 - dt2) < 1e-6) ? 0.5 * (dt1 + dt2) : (dt1 - dt2) / std::log(dt1 / dt2);

    // 3. Coeficientes de transferencia (Doran Fig. 9.17 / Prob. 9.13 + Gnielinski)
    double prCaldo = (p.cpCaldo * p.muCaldo) / p.kCaldo;
    double nuCaldo = 0.87 * std::pow(r.reAgitacion, 0.62) * std::pow(prCaldo, 1.0 / 3.0);
    r.ho = (nuCaldo * p.kCaldo) / p.diametroTanque;

    r.flujoRefReq = r.qTotal / (p.cpRef * (p.tRefOut - p.tRefIn));
    r.caudalSuficiente = (p.flujoRefDisp >= r.flujoRefReq);

    double areaTransversal = PI * std::pow(dIntTubo, 2) / 4.0;
    r.velRef = r.flujoRefReq / (p.rhoRef * areaTransversal);
    r.reRef = (p.rhoRef * r.velRef * dIntTubo) / p.muRef;
    double prRef = (p.cpRef * p.muRef) / p.kRef;

    r.diametroHelice = p.diametroTanque - 2.0 * p.holguraPared - p.dExtTubo;
    if (r.diametroHelice <= dIntTubo) {
        r.valido = false;
        r.mensajeError = "Error geometrico: La holgura radial es demasiado grande para el tanque.";
        return r;
    }

    double fDarcy = 0.0;
    double nuInt = 0.0;
    if (r.reRef < 2300.0) {
        fDarcy = 64.0 / r.reRef;
        nuInt = 3.66;
    } else {
        fDarcy = std::pow(0.79 * std::log(r.reRef) - 1.64, -2.0);
        double num = (fDarcy / 8.0) * (r.reRef - 1000.0) * prRef;
        double den = 1.0 + 12.7 * std::sqrt(fDarcy / 8.0) * (std::pow(prRef, 2.0 / 3.0) - 1.0);
        double nuRecto = num / den;
        nuInt = nuRecto * (1.0 + 3.5 * (dIntTubo / r.diametroHelice));
    }
    r.hi = (nuInt * p.kRef) / dIntTubo;

    double resTotal = (1.0 / r.ho) + p.factorEnsuciamiento +
                      ((p.dExtTubo * std::log(p.dExtTubo / dIntTubo)) / (2.0 * p.kPared)) +
                      (p.dExtTubo / (dIntTubo * r.hi));
    r.uGlobal = 1.0 / resTotal;

    // 4. Geometria del serpentin e hidraulica
    r.areaReq = r.qTotal / (r.uGlobal * r.dtlm);
    r.longitudTubo = r.areaReq / (PI * p.dExtTubo);

    double longPorVuelta = std::sqrt(std::pow(PI * r.diametroHelice, 2) + std::pow(p.pasoVueltas, 2));
    r.numVueltas = r.longitudTubo / longPorVuelta;
    r.alturaSerpentin = r.numVueltas * p.pasoVueltas;
    r.alturaLiquido = (4.0 * p.volumenLiquido) / (PI * std::pow(p.diametroTanque, 2));
    r.alturaTanqueEst = std::max(r.alturaLiquido * 1.25, p.diametroTanque * 1.5);

    double factorIto = 1.0 + 0.0823 * std::pow(dIntTubo / r.diametroHelice, 0.53) * std::pow(r.reRef, 0.25);
    r.caidaPresion = ((fDarcy * factorIto) * (r.longitudTubo / dIntTubo) * (0.5 * p.rhoRef * std::pow(r.velRef, 2))) / 1000.0;

    r.cabeEnLiquido = (r.alturaSerpentin < r.alturaLiquido);
    r.libreImpulsor = (r.diametroHelice > p.diametroImpulsor * 1.15);

    return r;
}

std::string escaparJson(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '"') out += "\\\"";
        else if (c == '\\') out += "\\\\";
        else out += c;
    }
    return out;
}

std::string resultadosAJson(const ParametrosProceso& p, const ResultadosDiseno& r) {
    std::ostringstream ss;
    ss.imbue(std::locale::classic());
    ss << std::fixed << std::setprecision(4);
    ss << "{"
       << "\"valido\":" << (r.valido ? "true" : "false") << ","
       << "\"mensajeError\":\"" << escaparJson(r.mensajeError) << "\","
       << "\"dhReaccion\":" << r.dhReaccion << ","
       << "\"qReaccionKw\":" << (r.qReaccion / 1000.0) << ","
       << "\"reAgitacion\":" << r.reAgitacion << ","
       << "\"potenciaNoAirKw\":" << (r.potenciaNoAir / 1000.0) << ","
       << "\"potenciaAireadaKw\":" << (r.potenciaAireada / 1000.0) << ","
       << "\"qTotalKw\":" << (r.qTotal / 1000.0) << ","
       << "\"dtlm\":" << r.dtlm << ","
       << "\"ho\":" << r.ho << ","
       << "\"flujoRefReq\":" << r.flujoRefReq << ","
       << "\"velRef\":" << r.velRef << ","
       << "\"reRef\":" << r.reRef << ","
       << "\"hi\":" << r.hi << ","
       << "\"uGlobal\":" << r.uGlobal << ","
       << "\"areaReq\":" << r.areaReq << ","
       << "\"longitudTubo\":" << r.longitudTubo << ","
       << "\"diametroHelice\":" << r.diametroHelice << ","
       << "\"numVueltas\":" << r.numVueltas << ","
       << "\"alturaSerpentin\":" << r.alturaSerpentin << ","
       << "\"alturaLiquido\":" << r.alturaLiquido << ","
       << "\"alturaTanqueEst\":" << r.alturaTanqueEst << ","
       << "\"diametroTanque\":" << p.diametroTanque << ","
       << "\"diametroImpulsor\":" << p.diametroImpulsor << ","
       << "\"dExtTubo\":" << p.dExtTubo << ","
       << "\"pasoVueltas\":" << p.pasoVueltas << ","
       << "\"caidaPresion\":" << r.caidaPresion << ","
       << "\"cabeEnLiquido\":" << (r.cabeEnLiquido ? "true" : "false") << ","
       << "\"libreImpulsor\":" << (r.libreImpulsor ? "true" : "false") << ","
       << "\"caudalSuficiente\":" << (r.caudalSuficiente ? "true" : "false")
       << "}";
    return ss.str();
}

ParametrosProceso parsearParametros(const std::string& target) {
    ParametrosProceso p;
    size_t qPos = target.find('?');
    if (qPos == std::string::npos) return p;

    std::string query = target.substr(qPos + 1);
    std::istringstream stream(query);
    std::string par;
    std::map<std::string, double> kv;

    while (std::getline(stream, par, '&')) {
        size_t eq = par.find('=');
        if (eq != std::string::npos) {
            std::string k = par.substr(0, eq);
            std::string v = par.substr(eq + 1);
            try { kv[k] = std::stod(v); } catch (...) {}
        }
    }

    auto get = [&](const std::string& key, double& dest) {
        if (kv.count(key)) dest = kv[key];
    };

    if (kv.count("modoEntalpia")) p.modoEntalpia = static_cast<int>(kv["modoEntalpia"]);
    get("dhCombSustrato", p.dhCombSustrato);
    get("dhCombBiomasa", p.dhCombBiomasa);
    get("rendimientoYxs", p.rendimientoYxs);
    get("dhDirecto", p.dhDirecto);
    get("tasaMolar", p.tasaMolar);
    get("factorCalor", p.factorCalor);

    get("diametroTanque", p.diametroTanque);
    get("volumenLiquido", p.volumenLiquido);
    get("diametroImpulsor", p.diametroImpulsor);
    get("velocidadAgitacion", p.velocidadAgitacion);
    get("numeroPotencia", p.numeroPotencia);
    get("flujoAire", p.flujoAire);
    get("rhoCaldo", p.rhoCaldo);
    get("muCaldo", p.muCaldo);
    get("kCaldo", p.kCaldo);
    get("cpCaldo", p.cpCaldo);

    get("tCaldo", p.tCaldo);
    get("tRefIn", p.tRefIn);
    get("tRefOut", p.tRefOut);
    get("flujoRefDisp", p.flujoRefDisp);
    get("cpRef", p.cpRef);
    get("rhoRef", p.rhoRef);
    get("muRef", p.muRef);
    get("kRef", p.kRef);

    get("dExtTubo", p.dExtTubo);
    get("espesorPared", p.espesorPared);
    get("kPared", p.kPared);
    get("factorEnsuciamiento", p.factorEnsuciamiento);
    get("holguraPared", p.holguraPared);
    get("pasoVueltas", p.pasoVueltas);

    return p;
}

const std::string PAGINA_PLANO = R"PLANO(
<!doctype html><html lang="es"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"><title>Plano en sección del biorreactor</title><style>
.iso{
  --anatomy-page:#0b0b0b;--anatomy-card:#131313;--anatomy-text:#a3a3a3;--anatomy-strong:#e6e6e6;--anatomy-muted:#6e6e6e;
  --anatomy-paper:#131313;--anatomy-face:#1b1b1b;--anatomy-top:#1c1c1c;
  --anatomy-shade-0:#0f0f0f;--anatomy-shade-1:#121212;--anatomy-shade-2:#151515;--anatomy-shade-3:#181818;
  --anatomy-lit-top:#262626;--anatomy-lit-shade:#1d1d1d;
  --anatomy-hi:#d4d4d4;--anatomy-mid:#6a6a6a;--anatomy-lo:#3e3e3e;--anatomy-faint:#2a2a2a;--anatomy-lit:#ffffff;
  --anatomy-halo:#050505;--anatomy-halo-opacity:.9;--anatomy-dot:#e6e6e6;
  --anatomy-red:#ff6b6b;--anatomy-green:#6be0a0;--anatomy-blue:#6b9bff;
  --anatomy-weight:.6px;--anatomy-ease:cubic-bezier(.32,.72,0,1);
  --anatomy-mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  --anatomy-lift:inset 0 1px 0 0 #ffffff0a,inset 0 0 0 1px #ffffff06,0 0 0 1px #00000026,0 1px 1px -.5px #00000014,0 3px 3px -1.5px #00000014,0 6px 6px -3px #00000014,0 12px 12px -6px #00000014;
}
.iso[data-theme="light"]{
  --anatomy-page:#f2f2ef;--anatomy-card:#fbfbfa;--anatomy-text:#57575a;--anatomy-strong:#18181a;--anatomy-muted:#8c8c90;
  --anatomy-paper:#fbfbfa;--anatomy-face:#ffffff;--anatomy-top:#ffffff;
  --anatomy-shade-0:#dcdcd9;--anatomy-shade-1:#e5e5e2;--anatomy-shade-2:#ededea;--anatomy-shade-3:#f4f4f1;
  --anatomy-lit-top:#ffffff;--anatomy-lit-shade:#eeeeeb;
  --anatomy-hi:#26262a;--anatomy-mid:#8a8a8f;--anatomy-lo:#b4b4b8;--anatomy-faint:#d4d4d6;--anatomy-lit:#000000;
  --anatomy-halo:#8a8a80;--anatomy-halo-opacity:.32;--anatomy-dot:#3a3a3e;
  --anatomy-red:#e5484d;--anatomy-green:#30a46c;--anatomy-blue:#3e63dd;
  --anatomy-lift:inset 0 1px 0 0 #ffffff,0 0 0 1px #0000000f,0 1px 1px -.5px #0000000a,0 3px 3px -1.5px #0000000a,0 6px 6px -3px #0000000a,0 12px 12px -6px #0000000a;
}
.iso-svg{display:block;width:100%;height:auto;overflow:visible}
.iso-svg path,.iso-svg circle{vector-effect:non-scaling-stroke}
.iso-fill{fill:var(--anatomy-paper);stroke:none}
.iso-solid[data-lit]>.iso-fill{fill:var(--anatomy-face)}
.iso-shade,.iso-top{stroke:none}
.iso-shade[data-shade="0"]{fill:var(--anatomy-shade-0)}
.iso-shade[data-shade="1"]{fill:var(--anatomy-shade-1)}
.iso-shade[data-shade="2"]{fill:var(--anatomy-shade-2)}
.iso-shade[data-shade="3"]{fill:var(--anatomy-shade-3)}
.iso-top{fill:var(--anatomy-top)}
.iso-solid[data-lit]>.iso-top{fill:var(--anatomy-lit-top)}
.iso-solid[data-lit]>.iso-shade{fill:var(--anatomy-lit-shade)}
.iso-solid[data-seamless]>.iso-shade,.iso-solid[data-seamless]>.iso-top{stroke-width:.4;stroke-linejoin:round}
.iso-solid[data-seamless]>.iso-shade[data-shade="0"]{stroke:var(--anatomy-shade-0)}
.iso-solid[data-seamless]>.iso-shade[data-shade="1"]{stroke:var(--anatomy-shade-1)}
.iso-solid[data-seamless]>.iso-shade[data-shade="2"]{stroke:var(--anatomy-shade-2)}
.iso-solid[data-seamless]>.iso-shade[data-shade="3"]{stroke:var(--anatomy-shade-3)}
.iso-solid[data-seamless]>.iso-top{stroke:var(--anatomy-top)}
.iso-solid[data-seamless][data-lit]>.iso-top{stroke:var(--anatomy-lit-top)}
.iso-solid[data-seamless][data-lit]>.iso-shade{stroke:var(--anatomy-lit-shade)}
.iso-solid[data-lit]>.iso-edge{stroke:var(--anatomy-lit)}
.iso-halo{fill:var(--anatomy-halo);stroke:none;filter:blur(9px);opacity:var(--anatomy-halo-opacity);transform:translateY(3px)}
.iso-line{fill:none;stroke:var(--anatomy-lo);stroke-width:var(--anatomy-weight);stroke-linecap:round;stroke-linejoin:round}
.iso-line[data-tone="hi"]{stroke:var(--anatomy-hi)}
.iso-line[data-tone="mid"]{stroke:var(--anatomy-mid)}
.iso-line[data-tone="lo"]{stroke:var(--anatomy-lo)}
.iso-line[data-tone="faint"]{stroke:var(--anatomy-faint)}
.iso-line[data-tone="lit"]{stroke:var(--anatomy-lit)}
.iso-line[data-tone="red"]{stroke:var(--anatomy-red)}
.iso-line[data-tone="green"]{stroke:var(--anatomy-green)}
.iso-line[data-tone="blue"]{stroke:var(--anatomy-blue)}
.iso-line[data-dotted]{stroke-dasharray:.01 3.5;stroke-width:1}
.iso-line[data-dashed]{stroke-dasharray:3 3}
.iso-line[data-flow]{stroke-dasharray:2 6}
.iso-wire-end{fill:var(--anatomy-paper);stroke:var(--anatomy-lo);stroke-width:var(--anatomy-weight)}
.iso-wire[data-tone="hi"] .iso-wire-end{stroke:var(--anatomy-hi)}
.iso-wire[data-tone="mid"] .iso-wire-end{stroke:var(--anatomy-mid)}
.iso-wire[data-tone="faint"] .iso-wire-end{stroke:var(--anatomy-faint)}
.iso-wire[data-tone="lit"] .iso-wire-end{stroke:var(--anatomy-lit)}
.iso-dots circle{fill:var(--anatomy-dot);stroke:none;opacity:.32}
.iso-dots[data-tone="hi"] circle{opacity:.85}
.iso-dots[data-tone="mid"] circle{opacity:.45}
.iso-dots[data-tone="lo"] circle{opacity:.22}
.iso-dots[data-tone="faint"] circle{opacity:.12}
.iso-face-text{fill:var(--anatomy-mid);font-family:var(--anatomy-mono);letter-spacing:.04em}
.iso-face-text[data-tone="hi"]{fill:var(--anatomy-hi)}
.iso-face-text[data-tone="lo"]{fill:var(--anatomy-lo)}
.iso-face-text[data-tone="faint"]{fill:var(--anatomy-faint)}
.tb-body{fill:var(--anatomy-shade-2);stroke:none}
.tb-shine{fill:var(--anatomy-top);stroke:none}
.tb-shade{fill:var(--anatomy-shade-0);stroke:none}
.tb-ring{fill:none;stroke:var(--anatomy-lo);stroke-width:.45;stroke-linecap:butt;stroke-linejoin:round}
.tb-ring[data-tone="mid"]{stroke:var(--anatomy-mid)}
.tb-ring[data-tone="faint"]{stroke:var(--anatomy-faint)}
.tb-ring[data-tone="hi"]{stroke:var(--anatomy-hi)}
@keyframes iso-flow{to{stroke-dashoffset:-32}}
@keyframes iso-pulse{0%,100%{opacity:.18}12%{opacity:1}40%{opacity:.18}}
.iso-line[data-flow]{animation:iso-flow 1600ms linear infinite}
.iso-dots circle[data-pulse]{animation:iso-pulse 1600ms steps(1,end) infinite;animation-delay:var(--pulse-at,0ms)}
.iso-figure{margin:32px 0 36px}
.iso-plate{position:relative;padding:52px 24px;border-radius:0;background:var(--anatomy-card);box-shadow:0 0 0 1px var(--anatomy-mid);--anatomy-paper:var(--anatomy-card)}
.iso-plate-corners{position:absolute;inset:18px 20px;pointer-events:none}
.iso-plate-corner{position:absolute;color:var(--anatomy-muted);font-family:var(--anatomy-mono);font-size:10.5px;line-height:1;letter-spacing:.04em;white-space:nowrap}
.iso-plate-corner[data-corner="fig"]{top:0;left:0;color:var(--anatomy-text)}
.iso-plate-corner[data-corner="title"]{top:0;right:0;text-transform:uppercase;letter-spacing:.1em}
.iso-plate-corner[data-corner="hint"]{bottom:0;left:0;text-transform:uppercase;letter-spacing:.1em}
.iso-plate-corner[data-corner="readout"]{right:0;bottom:0;color:var(--anatomy-text);font-variant-numeric:tabular-nums}
.iso-plate .iso-svg{margin:0 auto}
.iso-legend{display:grid;gap:8px;margin-top:12px;color:var(--anatomy-muted);font-size:13px;line-height:1.55}
.iso-keys{display:flex;flex-wrap:wrap;gap:6px 18px}
.iso-key{display:inline-flex;align-items:center;gap:7px}
.iso-swatch{width:16px;height:10px;flex:none;overflow:visible}
.iso-swatch path{fill:none;stroke:var(--anatomy-hi);stroke-width:1;vector-effect:non-scaling-stroke;stroke-linejoin:round;stroke-linecap:round}
.iso-swatch[data-mark="raised"] path,.iso-swatch[data-mark="flat"] path{stroke:var(--anatomy-mid)}
.iso-swatch[data-mark="lit"] path{stroke:var(--anatomy-lit)}
.iso-swatch[data-mark="dim"] path{stroke:var(--anatomy-lo)}
.iso-swatch[data-mark="dotted"] path{stroke-dasharray:.01 3;stroke-width:1.4}
.iso-swatch[data-mark="dashed"] path{stroke-dasharray:3 2.5}
.iso-swatch[data-mark="red"] path{stroke:var(--anatomy-red)}
.iso-swatch[data-mark="green"] path{stroke:var(--anatomy-green)}
.iso-swatch[data-mark="blue"] path{stroke:var(--anatomy-blue)}
.iso-swatch-crease{stroke:var(--anatomy-lo)!important}
@media (max-width:560px){.iso-plate{padding:48px 12px 70px}.iso-plate-corner[data-corner="hint"]{bottom:18px}}
@media (prefers-reduced-motion:reduce){.iso-line[data-flow],.iso-dots circle[data-pulse]{animation:none}}

*{box-sizing:border-box}
html{color-scheme:dark}body[data-theme="light"]{color-scheme:light}
body{margin:0;background:var(--anatomy-page);color:var(--anatomy-text);font:14px/1.5 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;-webkit-font-smoothing:antialiased;-webkit-tap-highlight-color:transparent}
.iso{--acc:#f5a524;--bx-edge:#0a0c0f;--bx-steel:#8793a1;--bx-shine:#e3eaf2;--bx-medium:#fb923c}
.iso[data-theme="light"]{--acc:#d97706;--bx-edge:#3b4148;--bx-steel:#c3cad3;--bx-shine:#ffffff}
button{font:inherit;color:inherit}
.app{max-width:1180px;margin:0 auto;padding:10px 12px 40px}
.top{display:flex;align-items:center;justify-content:space-between;gap:10px;padding:6px 2px 10px}
.top h1{font-size:16px;margin:0;color:var(--anatomy-strong);font-weight:600;letter-spacing:.01em}
.top p{margin:0;font-size:12px;color:var(--anatomy-muted)}
.btn{background:transparent;border:1px solid var(--anatomy-lo);border-radius:6px;padding:7px 11px;min-height:36px;cursor:pointer;font-size:12.5px}
.btn[aria-pressed="true"]{border-color:var(--acc);color:var(--anatomy-strong)}
.layout{display:grid;gap:14px}.layout>*{min-width:0}
@media(min-width:980px){.layout{grid-template-columns:minmax(0,1.05fr) minmax(0,.95fr);align-items:start}.stage{position:sticky;top:10px}}
.iso-plate.stagecard{padding:40px 6px 44px}
.iso-plate-corners{inset:14px 14px}
#fig{touch-action:manipulation}
#svg{max-height:82vh}
.bar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;justify-content:space-between;margin-top:10px}
.seg{display:inline-flex;border:1px solid var(--anatomy-lo);border-radius:6px;overflow:hidden}
.seg button{background:transparent;border:0;padding:7px 11px;min-height:36px;cursor:pointer;font-size:12.5px}
.seg button+button{border-left:1px solid var(--anatomy-lo)}
.seg button[aria-pressed="true"]{background:var(--anatomy-faint);color:var(--anatomy-strong)}
.tshow{margin-top:10px;display:grid;grid-template-columns:auto 1fr auto;gap:10px;align-items:center;font-size:12.5px}
.tshow output{font-variant-numeric:tabular-nums;color:var(--anatomy-strong);min-width:64px;text-align:right}
input[type=range]{width:100%;accent-color:var(--acc);min-height:28px}
.tabs{display:flex;border-bottom:1px solid var(--anatomy-lo);position:sticky;top:0;background:var(--anatomy-page);z-index:2}
.tab{flex:1;background:transparent;border:0;border-bottom:2px solid transparent;padding:11px 6px;cursor:pointer;font-size:13.5px;color:var(--anatomy-muted)}
.tab[aria-selected="true"]{color:var(--anatomy-strong);border-bottom-color:var(--acc)}
.panel{padding:14px 2px}
.chips{display:flex;gap:6px;overflow-x:auto;padding-bottom:8px;margin-bottom:10px;scrollbar-width:thin}
.chip{flex:none;background:transparent;border:1px solid var(--anatomy-lo);border-radius:999px;padding:6px 12px;min-height:34px;cursor:pointer;font-size:12.5px;white-space:nowrap}
.chip[aria-pressed="true"]{border-color:var(--acc);color:var(--anatomy-strong);background:var(--anatomy-faint)}
#info h3{margin:0 0 8px;font-size:17px;color:var(--anatomy-strong)}
#info h4,#calc h4{margin:16px 0 6px;font-size:12px;text-transform:uppercase;letter-spacing:.08em;color:var(--anatomy-muted);font-weight:600}
.facts{margin:0}.facts dt{font-size:12px;color:var(--anatomy-muted);margin-top:10px}.facts dd{margin:2px 0 0;color:var(--anatomy-text)}
.plist{list-style:none;margin:0;padding:0}.plist li{display:flex;justify-content:space-between;gap:10px;padding:7px 0;border-bottom:1px solid var(--anatomy-faint)}.plist b{color:var(--anatomy-strong);font-variant-numeric:tabular-nums;font-weight:600}
.hint,.muted{color:var(--anatomy-muted)}
fieldset{border:1px solid var(--anatomy-faint);border-radius:8px;margin:0 0 14px;padding:6px 12px 12px}
legend{padding:0 6px;font-size:12px;text-transform:uppercase;letter-spacing:.08em;color:var(--anatomy-muted)}
.field{padding:8px 0;border-bottom:1px solid var(--anatomy-faint)}.field:last-child{border-bottom:0}
.field label{display:block;font-size:13px;color:var(--anatomy-strong);margin-bottom:4px}
.inp{display:flex;align-items:center;gap:8px}
.inp input{flex:1;min-width:0;background:var(--anatomy-card);color:var(--anatomy-strong);border:1px solid var(--anatomy-lo);border-radius:6px;padding:8px 10px;min-height:40px;font:inherit;font-variant-numeric:tabular-nums}
.inp input:focus{outline:2px solid var(--acc);outline-offset:1px}
.u{font-size:12px;color:var(--anatomy-muted);min-width:64px}
.note{display:block;font-size:12px;color:var(--anatomy-muted);margin-top:2px}
.err{font-size:12.5px;margin-top:3px}.err:empty{display:none}.err[data-kind="bad"]{color:var(--anatomy-red)}.err[data-kind="warn"]{color:var(--acc)}
.field.bad .inp input{border-color:var(--anatomy-red)}.field.warn .inp input{border-color:var(--acc)}
#banner{border:1px solid var(--anatomy-red);color:var(--anatomy-red);border-radius:6px;padding:8px 10px;margin-bottom:12px;font-size:13px}
.ill{border-left:2px solid var(--acc);padding:4px 10px;margin:0 0 12px;font-size:12.5px;color:var(--anatomy-muted)}
.row{display:grid;grid-template-columns:1fr auto;gap:2px 10px;padding:8px 0;border-bottom:1px solid var(--anatomy-faint)}
.row .k{color:var(--anatomy-strong)}.row .v{text-align:right;font-variant-numeric:tabular-nums;color:var(--anatomy-strong);font-weight:600}
.row .f{grid-column:1/-1;font:11.5px var(--anatomy-mono);color:var(--anatomy-muted)}
.status{margin:10px 0 0;padding:8px 10px;border-radius:6px;border:1px solid var(--anatomy-lo)}
.status[data-s="coherente"]{border-color:var(--anatomy-green)}.status[data-s="revisar"]{border-color:var(--acc)}.status[data-s="inconsistente"]{border-color:var(--anatomy-red)}
details{margin-top:16px}summary{cursor:pointer;font-weight:600;color:var(--anatomy-strong)}
.sup{padding-left:18px;color:var(--anatomy-text)}.sup li{margin:5px 0}
.errs{padding-left:18px;color:var(--anatomy-red)}

.bp-grid path{stroke:var(--bp-grid);stroke-width:.5;fill:none}
.iso-svg .bp-o{fill:none;stroke:var(--bp-line);stroke-width:1.1;vector-effect:none}
.iso-svg .bp-o.bp-f{fill:var(--bp-fill)}
.iso-svg .bp-o.bp-cut{fill:url(#hatch)}
.bp-liq{fill:var(--bx-medium);opacity:.2;stroke:none;vector-effect:none}
.bp-hatch,.bp-thinline,.bp-ext,.bp-dimline{fill:none;stroke:var(--bp-dim);stroke-width:.7;vector-effect:none}
.bp-thinline{stroke:var(--bp-line);stroke-opacity:.6}
.bp-ext{stroke-opacity:.7}.bp-dimline{stroke-width:.9}
.bp-dimtext{fill:var(--bp-dim);font:600 11px var(--anatomy-mono);paint-order:stroke;stroke:var(--bp-paper);stroke-width:3px;stroke-linejoin:round}
.bp-tube,.bp-tube-edge,.bp-core,.bp-glow,.bp-flowline{fill:none;stroke-linecap:round;stroke-linejoin:round;vector-effect:none}
.bp-tube{stroke:var(--bp-line)}.bp-tube-edge{stroke:var(--bp-paper)}
.bp-glow{opacity:.12}.bp-core{opacity:.95}
.bp-flowline{stroke:#fff;stroke-opacity:.85;stroke-dasharray:2.4 9;animation:bx-flow var(--bx-flow-dur,2s) linear infinite}
.bp-hit{fill:none;stroke:transparent;pointer-events:stroke;stroke-linecap:round;vector-effect:none}
.bp-loop,.bp-heat,.bp-medline{fill:none;stroke-width:1.1;stroke-dasharray:3 5;stroke-linecap:round;animation:iso-flow 1.6s linear infinite}
.bp-loop{stroke:var(--bp-dim);stroke-opacity:.7}.bp-heat{stroke:var(--bx-medium);stroke-width:1.6;filter:saturate(1.6)}.bp-medline{stroke:var(--acc)}
.bp-heatg,.bp-dims{pointer-events:none}
.bp-blade{fill:var(--bp-fill);stroke:var(--bp-line);stroke-width:1;vector-effect:none}
.bp-dim-only{opacity:.25;pointer-events:none}
/* ilustración */
.bx-pipe path{vector-effect:none;fill:none;stroke-linecap:round;stroke-linejoin:round}
.bx-edge{stroke:var(--bx-edge)}.bx-body{stroke:var(--bx-steel)}.bx-shine{stroke:var(--bx-shine);opacity:.75}
.bx-core{opacity:.92}.bx-glow{stroke-linecap:round}
.bx-flowline{stroke:#fff;stroke-opacity:.8;stroke-dasharray:2.4 9;animation:bx-flow var(--bx-flow-dur,2s) linear infinite}
@keyframes bx-flow{to{stroke-dashoffset:-11.4}}
.bx-hit{fill:none;stroke:transparent;pointer-events:stroke;stroke-linecap:round;vector-effect:none}
[data-part]{cursor:pointer}
[data-part].sel{filter:drop-shadow(0 0 2.5px var(--acc)) drop-shadow(0 0 .8px var(--acc))}
.bx-liq{fill:var(--bx-medium);stroke:none;opacity:.15;pointer-events:all}
.bx-liq-front{opacity:.07;pointer-events:none}.bx-liq-surf{opacity:.18;pointer-events:none}
.bx-veil{pointer-events:none}
.bx-swirl{fill:none;stroke:var(--anatomy-hi);stroke-opacity:.4;stroke-width:.8;animation:iso-flow 3s linear infinite}
.bx-loop{stroke:var(--anatomy-mid)}
.bx-glass{opacity:.2;pointer-events:none}
.bx-heat{pointer-events:none}
.bx-heat-line{fill:none;stroke:var(--bx-medium);stroke-width:1.5;stroke-dasharray:2 4;stroke-linecap:round;animation:iso-flow 1.5s linear infinite;filter:saturate(1.6)}
.bx-heat-head{fill:none;stroke:var(--bx-medium);stroke-width:1.2;stroke-linecap:round;stroke-linejoin:round;filter:saturate(1.6)}
.bx-med-line,.bx-cool-line{fill:none;stroke-width:1.1;stroke-dasharray:2 5;stroke-linecap:round;animation:iso-flow 1.6s linear infinite}
.bx-med-line{stroke:var(--acc)}
.bx-label{cursor:pointer}
.bx-label text{fill:var(--anatomy-text);font-family:var(--anatomy-mono);paint-order:stroke;stroke:var(--anatomy-card);stroke-width:3px;stroke-linejoin:round;letter-spacing:.02em}
.bx-leader{fill:none;stroke:var(--anatomy-mid);stroke-width:.6}.bx-leader-dot{fill:var(--anatomy-hi)}
.bx-label.sel text{fill:var(--acc);font-weight:700}.bx-label.sel .bx-leader{stroke:var(--acc)}
.paused *{animation-play-state:paused!important}
@media(prefers-reduced-motion:reduce){.bx-flowline,.bx-swirl,.bx-heat-line,.bx-med-line,.bx-cool-line{animation:none}}
body{--bp-paper:#0c2a4d;--bp-line:#d8ecff;--bp-dim:#7fd0ff;--bp-grid:#16406f;--bp-fill:#12365f;--acc:#ffd166}body[data-theme="light"]{--bp-paper:#fbfbf8;--bp-line:#1d2b3a;--bp-dim:#0b6bcb;--bp-grid:#e6e9ee;--bp-fill:#e9eef5;--acc:#d9480f}.iso-plate.stagecard{background:var(--bp-paper);border-color:var(--bp-dim)}.iso-plate-corner{color:var(--bp-dim)!important}</style><link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,300..600;1,9..144,300..500&family=JetBrains+Mono:wght@400;500;600&family=Spline+Sans:wght@400;500;600&display=swap" rel="stylesheet">
<style id="amber-plano-theme">
  :root, .iso {
    --ink: #16201f;
    --abyss: #0e2a29;
    --teal: #1f6b6a;
    --teal-bright: #2d8b8b;
    --seafoam: #9ec9c6;
    --parchment: #f4f1e8;
    --parchment-deep: #ebe6d6;
    --cream: #fbfaf4;
    --gold: #b08440;
    --gold-soft: #c89a52;
    --rule: #d8d2bf;
    --code-bg: #10302e;
    --code-ink: #dfeae7;
    --anatomy-mono: 'JetBrains Mono', ui-monospace, "SF Mono", Menlo, Consolas, monospace;
    --acc: #e6b566;
  }
  /* Modo claro editorial (pergamino y tinta abisal) */
  .iso[data-theme="light"], body[data-theme="light"] {
    --anatomy-page: #f4f1e8;
    --anatomy-card: #fbfaf4;
    --anatomy-text: #3c4a48;
    --anatomy-strong: #0e2a29;
    --anatomy-muted: #5a6b69;
    --anatomy-paper: #fbfaf4;
    --anatomy-face: #ffffff;
    --anatomy-top: #ffffff;
    --anatomy-hi: #0e2a29;
    --anatomy-mid: #5a6b69;
    --anatomy-lo: #b8b19c;
    --anatomy-faint: #e5dfd0;
    --anatomy-lit: #1f6b6a;
    --anatomy-red: #b84a39;
    --anatomy-green: #1f6b6a;
    --anatomy-blue: #2d8b8b;
    --bp-paper: #fbfaf4;
    --bp-line: #0e2a29;
    --bp-dim: #1f6b6a;
    --bp-grid: #ebe6d6;
    --bp-fill: #f0ece1;
    --acc: #b08440;
  }
  /* Modo oscuro abisal (inspirado en bloques de codigo AMBER) */
  .iso[data-theme="dark"], body[data-theme="dark"] {
    --anatomy-page: #0b2221;
    --anatomy-card: #10302e;
    --anatomy-text: #bcd9d4;
    --anatomy-strong: #f4f1e8;
    --anatomy-muted: #7fa9a4;
    --anatomy-paper: #10302e;
    --anatomy-hi: #dfeae7;
    --anatomy-mid: #7fa9a4;
    --anatomy-lo: #2b5451;
    --anatomy-faint: #19423f;
    --bp-paper: #10302e;
    --bp-line: #dfeae7;
    --bp-dim: #9ec9c6;
    --bp-grid: #163d3a;
    --bp-fill: #153b38;
    --acc: #e6b566;
  }
  body {
    font-family: 'Spline Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif !important;
    background-color: var(--anatomy-page) !important;
    background-image:
      radial-gradient(circle at 12% 8%, rgba(45, 139, 139, 0.07), transparent 40%),
      radial-gradient(circle at 88% 92%, rgba(176, 132, 64, 0.07), transparent 42%);
  }
  .top h1, #info h3, .col-card h4 {
    font-family: 'Fraunces', Georgia, serif !important;
    font-weight: 500 !important;
    letter-spacing: -0.01em;
  }
  .top h1 { font-size: 1.35rem !important; }
  .iso-plate.stagecard {
    border-radius: 12px !important;
    box-shadow: 0 14px 32px -20px rgba(14, 42, 41, 0.35), 0 0 0 1px var(--bp-dim) !important;
  }
  fieldset, .status {
    background: var(--anatomy-card);
    border-color: var(--anatomy-lo) !important;
    border-radius: 10px !important;
  }
  legend, #info h4, #calc h4, .iso-plate-corner {
    font-family: 'JetBrains Mono', monospace !important;
    letter-spacing: 0.14em !important;
  }
  .chip, .btn, .seg {
    border-radius: 8px !important;
    transition: all 0.2s ease;
  }
  .chip { border-radius: 999px !important; }
</style>
</head>
<body class="iso" data-theme="dark">
<div id="nucoil-nav" style="background:rgba(244,241,232,0.92);backdrop-filter:blur(10px);border-bottom:1px solid #d8d2bf;padding:10px 18px;display:flex;justify-content:space-between;align-items:center;font-family:'Spline Sans',system-ui,sans-serif;font-size:13px;color:#0e2a29;position:sticky;top:0;z-index:40">
  <span style="font-family:'Fraunces',serif;font-size:15px;font-weight:500">NuCoil <em style="color:#1f6b6a;font-style:italic">· Plano en Sección</em></span>
  <a id="btn-volver-calc" href="/" style="background:#1f6b6a;color:#fbfaf4;text-decoration:none;padding:6px 14px;border-radius:8px;font-family:'JetBrains Mono',monospace;font-size:11.5px;letter-spacing:.04em">← Calculadora Térmica</a>
</div>
<div class="app">
<header class="top"><div><h1>Plano en sección del biorreactor</h1><p>Alzado con cotas que siguen los parámetros</p></div><button class="btn" id="btnTheme" type="button">Claro / oscuro</button></header>
<div class="layout">
<section class="stage"><figure class="iso-figure" style="margin:0"><div class="iso-plate stagecard">
<div class="iso-plate-corners" aria-hidden="true"><span class="iso-plate-corner" data-corner="fig">Lámina 4</span><span class="iso-plate-corner" data-corner="title">Sección A-A</span><span class="iso-plate-corner" data-corner="hint">Toca un componente</span><span class="iso-plate-corner" data-corner="readout" data-readout></span></div>
<div id="fig"><svg id="svg" class="iso-svg" viewBox="0 0 620 470" role="img" aria-label="Biorreactor cilíndrico vertical de acero inoxidable en corte parcial, con agitador de dos turbinas, cuatro deflectores y un serpentín helicoidal de enfriamiento. El refrigerante entra por abajo, sube por la hélice y sale por la tapa."><defs><pattern id="hatch" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="5" height="5" fill="var(--bp-fill)"/><path d="M0 0V5" stroke="var(--bp-line)" stroke-width="1" stroke-opacity=".8"/></pattern></defs><g id="g-base"></g><g id="g-back"></g><g id="g-coilBack"></g><g id="g-flow"></g><g id="g-dynBack"></g><g id="g-shaft"></g><g id="g-dyn"></g><g id="g-dynFront"></g><g id="g-front"></g><g id="g-heat"></g><g id="g-veil"></g><g id="g-shell"></g><g id="g-top"></g><g id="g-labels"></g></svg></div>
</div>
<figcaption class="iso-legend"><span class="iso-keys"><span class="iso-key"><svg class="iso-swatch" viewBox="0 0 16 10"><path d="M1.5 5h13" style="stroke:var(--anatomy-blue);stroke-width:2"/></svg>Refrigerante frío</span><span class="iso-key"><svg class="iso-swatch" viewBox="0 0 16 10"><path d="M1.5 5h13" style="stroke:var(--anatomy-red);stroke-width:2"/></svg>Medio caliente</span><span class="iso-key"><svg class="iso-swatch" viewBox="0 0 16 10"><path d="M1.5 5h13" style="stroke-dasharray:3 2.5"/></svg>Calor hacia el serpentín</span></span></figcaption></figure>
<div class="bar"><div class="seg" role="group" aria-label="Vista del tanque"><button type="button" data-view="cut" aria-pressed="true">Corte</button><button type="button" data-view="glass" aria-pressed="false">Transparente</button><button type="button" data-view="closed" aria-pressed="false">Cerrado</button></div><div style="display:flex;gap:8px"><button class="btn" id="btnLabels" type="button" aria-pressed="true">Etiquetas</button><button class="btn" id="btnPause" type="button" aria-pressed="false">Pausar</button></div></div>
<div class="tshow"><label for="tshow">Temperatura del cultivo mostrada</label><input id="tshow" type="range"><output id="tshowv"></output></div>
</section>
<section class="side"><div class="tabs" role="tablist"><button class="tab" role="tab" data-tab="comp" aria-selected="true">Componente</button><button class="tab" role="tab" data-tab="par" aria-selected="false">Parámetros</button><button class="tab" role="tab" data-tab="calc" aria-selected="false">Cálculos</button></div>
<div class="panel" data-tab="comp"><div class="chips" id="chips"></div><div id="info"></div></div>
<div class="panel" data-tab="par" hidden><p class="ill">Los valores iniciales son ilustrativos para poder dibujar el modelo. No son resultados experimentales.</p><div id="banner" hidden></div><div id="form"></div><button class="btn" id="btnReset" type="button">Restablecer valores ilustrativos</button></div>
<div class="panel" data-tab="calc" hidden><div id="calc"></div></div>
</section></div></div>
<script type="module">const AZIMUTH = 45;
const ELEVATION = 30;
const SHADES = 4;
const RADIANS = Math.PI / 180;
const UNIT = Math.sqrt(1.6);
function cameraOf({ scale, azimuth = AZIMUTH, elevation = ELEVATION }) {
    const a = azimuth * RADIANS;
    const e = elevation * RADIANS;
    return { sinA: Math.sin(a), cosA: Math.cos(a), sinE: Math.sin(e), cosE: Math.cos(e), k: scale * UNIT };
}
function iso([x, y, z], projection) {
    const { origin } = projection;
    const { sinA, cosA, sinE, cosE, k } = cameraOf(projection);
    return [origin[0] + (x * sinA - y * cosA) * k, origin[1] + ((x * cosA + y * sinA) * sinE - z * cosE) * k];
}
function towardViewer(projection) {
    const { sinA, cosA } = cameraOf(projection);
    return [cosA, sinA];
}
function depthOf([x, y, z], projection) {
    const { sinA, cosA, sinE, cosE } = cameraOf(projection);
    return (x * cosA + y * sinA) * cosE + z * sinE;
}
const tenth = (value) => {
    const text = (Math.round(value * 10) / 10).toFixed(1).replace(/\.0$/, "");
    return text === "-0" ? "0" : text;
};
function pathOf(points, closed = false) {
    if (!points.length)
        return "";
    return points.map(([x, y], at) => `${at ? "L" : "M"}${tenth(x)} ${tenth(y)}`).join("") + (closed ? "Z" : "");
}
const joinPaths = (paths) => paths.filter(Boolean).join("");
function roundedPlan({ x, y, w, d, r }, steps = 8) {
    const radius = Math.max(0, Math.min(r, w / 2, d / 2));
    const corners = [
        [x + w - radius, y + radius, -90],
        [x + w - radius, y + d - radius, 0],
        [x + radius, y + d - radius, 90],
        [x + radius, y + radius, 180],
    ];
    const out = [];
    for (const [cx, cy, start] of corners) {
        const count = radius > 0 ? steps : 0;
        for (let step = 0; step <= count; step++) {
            const angle = (start + (count ? (90 * step) / count : 45)) * RADIANS;
            out.push([cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)]);
        }
    }
    return out;
}
function signedArea(points) {
    let area = 0;
    for (let index = 0; index < points.length; index++) {
        const [ax, ay] = points[index];
        const [bx, by] = points[(index + 1) % points.length];
        area += ax * by - bx * ay;
    }
    return area / 2;
}
const oriented = (points) => (signedArea(points) < 0 ? [...points].reverse() : points);
function convexHull(points) {
    const sorted = [...points].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    if (sorted.length < 3)
        return sorted;
    const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const lower = [];
    for (const point of sorted) {
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], point) <= 0)
            lower.pop();
        lower.push(point);
    }
    const upper = [];
    for (let index = sorted.length - 1; index >= 0; index--) {
        const point = sorted[index];
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], point) <= 0)
            upper.pop();
        upper.push(point);
    }
    return [...lower.slice(0, -1), ...upper.slice(0, -1)];
}
function runsOf(flags, want) {
    const count = flags.length;
    const runs = [];
    const start = flags.findIndex((flag, index) => flag !== flags[(index + count - 1) % count]);
    if (start < 0)
        return flags[0] === want ? [Array.from({ length: count + 1 }, (_, index) => index % count)] : [];
    let run = null;
    for (let step = 0; step < count; step++) {
        const index = (start + step) % count;
        if (flags[index] === want) {
            if (!run)
                run = [index];
            run.push((index + 1) % count);
        }
        else if (run) {
            runs.push(run);
            run = null;
        }
    }
    if (run)
        runs.push(run);
    return runs;
}
function extrude(ring, z, height, projection, options = {}) {
    const convex = options.convex ?? true;
    const count = ring.length;
    const counter = signedArea(ring) > 0;
    const [vx, vy] = towardViewer(projection);
    const { sinA, cosA } = cameraOf(projection);
    const normals = ring.map(([ax, ay], index) => {
        const [bx, by] = ring[(index + 1) % count];
        const dx = bx - ax;
        const dy = by - ay;
        const length = Math.hypot(dx, dy) || 1;
        return (counter ? [dy / length, -dx / length] : [-dy / length, dx / length]);
    });
    const facing = normals.map(([nx, ny]) => nx * vx + ny * vy > 1e-9);
    const top = ring.map(([x, y]) => iso([x, y, z + height], projection));
    const bottom = ring.map(([x, y]) => iso([x, y, z], projection));
    const outline = [];
    const crease = [];
    for (const run of runsOf(facing, false))
        outline.push(pathOf(run.map((index) => top[index])));
    for (const run of runsOf(facing, true)) {
        crease.push(pathOf(run.map((index) => top[index])));
        if (height > 0)
            outline.push(pathOf(run.map((index) => bottom[index])));
    }
    if (height > 0) {
        for (let index = 0; index < count; index++) {
            if (facing[(index + count - 1) % count] !== facing[index])
                outline.push(pathOf([top[index], bottom[index]]));
        }
    }
    let fill;
    if (convex) {
        fill = pathOf(convexHull([...top, ...bottom]), true);
    }
    else {
        const parts = [pathOf(oriented(top), true), pathOf(oriented(bottom), true)];
        for (let index = 0; index < count; index++) {
            if (!facing[index])
                continue;
            const next = (index + 1) % count;
            parts.push(pathOf(oriented([top[index], top[next], bottom[next], bottom[index]]), true));
        }
        fill = parts.join("");
    }
    const buckets = Array.from({ length: SHADES }, () => []);
    if (height > 0) {
        const [lx, ly] = [-sinA, cosA];
        for (let index = 0; index < count; index++) {
            if (!facing[index])
                continue;
            const [nx, ny] = normals[index];
            const lit = (nx * lx + ny * ly + 1) / 2;
            const shade = Math.min(SHADES - 1, Math.max(0, Math.floor(lit * SHADES)));
            const next = (index + 1) % count;
            buckets[shade].push(pathOf(oriented([top[index], top[next], bottom[next], bottom[index]]), true));
        }
    }
    let bevel;
    if (options.bevel && options.bevel > 0) {
        const by = options.bevel;
        const inset = ring.map(([x, y], index) => {
            const [ax, ay] = normals[(index + count - 1) % count];
            const [bx, by2] = normals[index];
            const mx = ax + bx;
            const my = ay + by2;
            const length = Math.hypot(mx, my) || 1;
            return [x - (mx / length) * by, y - (my / length) * by];
        });
        const insetTop = inset.map(([x, y]) => iso([x, y, z + height], projection));
        bevel = runsOf(facing, true)
            .map((run) => pathOf([top[run[0]], ...run.slice(1, -1).map((index) => insetTop[index]), top[run[run.length - 1]]]))
            .join("");
    }
    return {
        fill,
        outline: outline.join(""),
        crease: crease.join(""),
        top: pathOf(oriented(top), true),
        shades: buckets.map((paths) => paths.join("")),
        bevel,
    };
}
const ringOf = (plan, steps = 8) => roundedPlan(plan, steps);
function slabOf(plan, z, height, projection, steps = 8, bevel = 0) {
    return extrude(ringOf(plan, steps), z, height, projection, { convex: true, bevel });
}
function circleRing(cx, cy, radius, steps = 40) {
    return Array.from({ length: steps }, (_, step) => {
        const angle = (step / steps) * Math.PI * 2;
        return [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)];
    });
}
function cylinder(cx, cy, radius, z, height, projection, steps = 32, bevel = 0) {
    return extrude(circleRing(cx, cy, radius, steps), z, height, projection, { convex: true, bevel });
}
const insetPlan = ({ x, y, w, d, r }, by) => ({ x: x + by, y: y + by, w: w - by * 2, d: d - by * 2, r: Math.max(0, r - by) });
function onTop(points, z, projection, closed = false) {
    return pathOf(points.map(([x, y]) => iso([x, y, z], projection)), closed);
}
const planOutline = (plan, z, projection, steps = 8) => onTop(ringOf(plan, steps), z, projection, true);
const lineOnTop = (from, to, z, projection) => onTop([from, to], z, projection);
const segment = (from, to, projection) => pathOf([iso(from, projection), iso(to, projection)]);
const rise = (at, from, to, projection) => segment([at[0], at[1], from], [at[0], at[1], to], projection);
const ring = (cx, cy, radius, z, projection, steps = 40) => onTop(circleRing(cx, cy, radius, steps), z, projection, true);
const pointsAt = (points, projection) => points.map((point) => iso(point, projection));
function dotGrid(x, y, columns, rows, gap, z, projection) {
    const dots = [];
    for (let row = 0; row < rows; row++) {
        for (let column = 0; column < columns; column++)
            dots.push(iso([x + column * gap, y + row * gap, z], projection));
    }
    return dots;
}
function corners({ x, y, w, d }, inset = 0) {
    return [
        [x + inset, y + inset],
        [x + w - inset, y + inset],
        [x + inset, y + d - inset],
        [x + w - inset, y + d - inset],
    ];
}
function boxCorners({ x, y, w, d }, z0, z1) {
    const out = [];
    for (const px of [x, x + w])
        for (const py of [y, y + d])
            for (const z of [z0, z1])
                out.push([px, py, z]);
    return out;
}
function knurl(cx, cy, radius, z0, z1, count, projection, inset = 0.4) {
    const [vx, vy] = towardViewer(projection);
    const lines = [];
    for (let index = 0; index < count; index++) {
        const angle = (index / count) * Math.PI * 2;
        const nx = Math.cos(angle);
        const ny = Math.sin(angle);
        if (nx * vx + ny * vy < 0.15)
            continue;
        lines.push(segment([cx + radius * nx, cy + radius * ny, z0 + inset], [cx + radius * nx, cy + radius * ny, z1 - inset], projection));
    }
    return lines.join("");
}
function sideArc(cx, cy, radius, z, projection, steps = 24) {
    const [vx, vy] = towardViewer(projection);
    const facing = Math.atan2(vy, vx);
    const points = [];
    for (let step = 0; step <= steps; step++) {
        const angle = facing - Math.PI / 2 + (step / steps) * Math.PI;
        points.push(iso([cx + radius * Math.cos(angle), cy + radius * Math.sin(angle), z], projection));
    }
    return pathOf(points);
}
function ringSeam(ring, z, projection) {
    const count = ring.length;
    const counter = signedArea(ring) > 0;
    const [vx, vy] = towardViewer(projection);
    const facing = ring.map(([ax, ay], index) => {
        const [bx, by] = ring[(index + 1) % count];
        const dx = bx - ax;
        const dy = by - ay;
        const [nx, ny] = counter ? [dy, -dx] : [-dy, dx];
        return nx * vx + ny * vy > 1e-9;
    });
    return runsOf(facing, true).map((run) => pathOf(run.map((index) => iso([ring[index][0], ring[index][1], z], projection)))).join("");
}
const sideSeam = (plan, z, projection, steps = 8) => ringSeam(ringOf(plan, steps), z, projection);
function sideRing(center, radius, axis, projection, steps = 20) {
    const [px, py, pz] = center;
    const points = [];
    for (let step = 0; step < steps; step++) {
        const angle = (step / steps) * Math.PI * 2;
        const u = radius * Math.cos(angle);
        const v = radius * Math.sin(angle);
        points.push(iso(axis === "x" ? [px, py + u, pz + v] : [px + u, py, pz + v], projection));
    }
    return pathOf(points, true);
}
function coil(from, to, center, radius, turns, axis, projection, perTurn = 16) {
    const count = Math.max(1, Math.round(turns * perTurn));
    const back = [];
    const front = [];
    let run = [];
    let near = null;
    for (let index = 0; index <= count; index++) {
        const share = index / count;
        const angle = share * turns * Math.PI * 2;
        const u = radius * Math.cos(angle);
        const v = radius * Math.sin(angle);
        const along = from + (to - from) * share;
        const point = axis === "x" ? [along, center[0] + u, center[1] + v] : [center[0] + u, along, center[1] + v];
        const normal = axis === "x" ? [0, Math.cos(angle), Math.sin(angle)] : [Math.cos(angle), 0, Math.sin(angle)];
        const facing = depthOf(normal, projection) > 0;
        const flat = iso(point, projection);
        if (near === null)
            near = facing;
        if (facing !== near) {
            run.push(flat);
            (near ? front : back).push(pathOf(run));
            run = [flat];
            near = facing;
            continue;
        }
        run.push(flat);
    }
    if (run.length > 1)
        (near ? front : back).push(pathOf(run));
    return { back: back.join(""), front: front.join("") };
}
function topTicks(from, to, step, every, edge, z, lengths, projection, axis = "x") {
    const minor = [];
    const major = [];
    const first = Math.ceil(from / step) * step;
    for (let at = first, index = Math.round(first / step); at <= to + 1e-9; at += step, index++) {
        const big = index % every === 0;
        const length = big ? lengths[1] : lengths[0];
        const a = axis === "x" ? [at, edge] : [edge, at];
        const b = axis === "x" ? [at, edge + length] : [edge + length, at];
        (big ? major : minor).push(lineOnTop(a, b, z, projection));
    }
    return { minor: minor.join(""), major: major.join("") };
}
function sideTicks(from, to, step, every, face, top, lengths, projection, axis = "x") {
    const minor = [];
    const major = [];
    const first = Math.ceil(from / step) * step;
    for (let at = first, index = Math.round(first / step); at <= to + 1e-9; at += step, index++) {
        const big = index % every === 0;
        const length = big ? lengths[1] : lengths[0];
        const a = axis === "x" ? [at, face, top] : [face, at, top];
        const b = axis === "x" ? [at, face, top - length] : [face, at, top - length];
        (big ? major : minor).push(segment(a, b, projection));
    }
    return { minor: minor.join(""), major: major.join("") };
}
function radialTicks(cx, cy, radius, count, every, z, lengths, projection, from = 0, span = Math.PI * 2) {
    const minor = [];
    const major = [];
    const closed = Math.abs(span - Math.PI * 2) < 1e-9;
    const steps = closed ? count : count + 1;
    for (let index = 0; index < steps; index++) {
        const angle = from + (index / count) * span;
        const big = index % every === 0;
        const inner = radius - (big ? lengths[1] : lengths[0]);
        const sx = Math.sin(angle);
        const sy = -Math.cos(angle);
        (big ? major : minor).push(lineOnTop([cx + sx * inner, cy + sy * inner], [cx + sx * radius, cy + sy * radius], z, projection));
    }
    return { minor: minor.join(""), major: major.join("") };
}
const haloOf = (plan, z, height, projection) => slabOf(plan, z, height, projection).fill;
function axisVector(projection, [dx, dy, dz]) {
    const origin = iso([0, 0, 0], projection);
    const moved = iso([dx, dy, dz], projection);
    return [moved[0] - origin[0], moved[1] - origin[1]];
}
const fine = (value) => Math.round(value * 100) / 100;
function translateAlong(projection, by) {
    const [x, y] = axisVector(projection, by);
    return `translate(${fine(x)} ${fine(y)})`;
}
const thousandth = (value) => Math.round(value * 1000) / 1000;
function axesAt(at, z, projection) {
    const origin = iso([at[0], at[1], z], projection);
    const along = (point) => {
        const [px, py] = iso(point, projection);
        return [px - origin[0], py - origin[1]];
    };
    return { origin, x: along([at[0] + 1, at[1], z]), y: along([at[0], at[1] + 1, z]), z: along([at[0], at[1], z + 1]) };
}
const matrixOf = (u, v, origin) => `matrix(${thousandth(u[0])} ${thousandth(u[1])} ${thousandth(v[0])} ${thousandth(v[1])} ${thousandth(origin[0])} ${thousandth(origin[1])})`;
function holdStill(element, [dx, dy]) {
    const text = dx || dy ? `translate(${-dx} ${-dy})` : "";
    element.querySelectorAll("[data-hold]").forEach((node) => (text ? node.setAttribute("transform", text) : node.removeAttribute("transform")));
}
function topMatrix(at, z, projection, along = "x") {
    const { origin, x, y } = axesAt(at, z, projection);
    return along === "x" ? matrixOf(x, y, origin) : matrixOf([-y[0], -y[1]], x, origin);
}
function sideMatrix(at, z, projection, face) {
    const { origin, x, y, z: up } = axesAt(at, z, projection);
    const down = [-up[0], -up[1]];
    return face === "left" ? matrixOf(x, down, origin) : matrixOf([-y[0], -y[1]], down, origin);
}
function screenRise(height, projection) {
    const { cosE, k } = cameraOf(projection);
    return height * cosE * k;
}
function fitProjection(points, width, height, { pad = 24, headroom = 0, azimuth = AZIMUTH, elevation = ELEVATION } = {}) {
    const unit = { origin: [0, 0], scale: 1, azimuth, elevation };
    const flat = points.map((point) => iso(point, unit));
    const xs = flat.map(([x]) => x);
    const ys = flat.map(([, y]) => y);
    const left = Math.min(...xs);
    const right = Math.max(...xs);
    const up = Math.min(...ys);
    const down = Math.max(...ys);
    const scale = Math.min((width - pad * 2) / (right - left), (height - pad * 2 - headroom) / (down - up));
    return {
        origin: [width / 2 - ((left + right) / 2) * scale, headroom + (height - headroom) / 2 - ((up + down) / 2) * scale],
        scale,
        azimuth,
        elevation,
    };
}
const escape = (text) => text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
const attrs = (pairs) => Object.entries(pairs)
    .filter(([, value]) => value !== undefined && value !== false)
    .map(([key, value]) => ` ${key}="${escape(String(value))}"`)
    .join("");
function solidSvg(paths, { tone = "hi", crease = "faint", halo = false, lit = false, flat = false, className, inner = "" } = {}) {
    const parts = [];
    if (halo)
        parts.push(`<path class="iso-halo" d="${paths.fill}"/>`);
    parts.push(`<path class="iso-fill" d="${paths.fill}"/>`);
    if (!flat) {
        paths.shades.forEach((d, shade) => {
            if (d)
                parts.push(`<path class="iso-shade" data-shade="${shade}" d="${d}"/>`);
        });
        parts.push(`<path class="iso-top" d="${paths.top}"/>`);
    }
    if (crease !== "none" && paths.crease)
        parts.push(`<path class="iso-line iso-crease" data-tone="${crease}" d="${paths.crease}"/>`);
    if (crease !== "none" && paths.bevel)
        parts.push(`<path class="iso-line iso-bevel" data-tone="${crease}" d="${paths.bevel}"/>`);
    if (inner)
        parts.push(inner);
    parts.push(`<path class="iso-line iso-edge" data-tone="${tone}" d="${paths.outline}"/>`);
    return `<g${attrs({ class: className ? `iso-solid ${className}` : "iso-solid", "data-lit": lit ? "" : undefined, "data-seamless": paths.seamless ? "" : undefined })}>${parts.join("")}</g>`;
}
function lineSvg(d, { tone = "lo", dotted = false, dashed = false, flow = false, className, free } = {}) {
    if (!d)
        return "";
    return `<path${attrs({
        class: className ? `iso-line ${className}` : "iso-line",
        "data-tone": tone,
        "data-free": free === true ? "" : free || undefined,
        "data-dotted": dotted ? "" : undefined,
        "data-dashed": dashed ? "" : undefined,
        "data-flow": flow ? "" : undefined,
        d,
    })}/>`;
}
function dotsSvg(points, { size = 0.5, tone = "mid", pulse = false } = {}) {
    const dots = points
        .map(([x, y, alpha], index) => `<circle cx="${tenth(x)}" cy="${tenth(y)}" r="${size}"${alpha !== undefined && alpha < 1 ? ` fill-opacity="${Math.max(0, alpha).toFixed(2)}"` : ""}${pulse ? ` data-pulse="" style="--pulse-at:${(index * 137) % 1600}ms"` : ""}/>`)
        .join("");
    return `<g class="iso-dots" data-tone="${tone}">${dots}</g>`;
}
function fadedSvg(list, style = {}) {
    const groups = new Map();
    for (const entry of list) {
        const [d, alpha] = Array.isArray(entry) ? entry : [entry.d, entry.alpha];
        if (!d || !(alpha > 0.02))
            continue;
        const level = Math.min(1, Math.round(alpha * 20) / 20);
        groups.set(level, (groups.get(level) ?? "") + d);
    }
    return [...groups].map(([alpha, d]) => (alpha >= 1 ? lineSvg(d, style) : `<g opacity="${alpha.toFixed(2)}">${lineSvg(d, style)}</g>`)).join("");
}
function fadeLineSvg(points, { fade = [0, 0], tone = "lo", className, stops = 6 } = {}) {
    if (points.length < 2)
        return "";
    const ease = (t) => t * t * (3 - 2 * t);
    const hashOf = (text) => {
        let hash = 2166136261;
        for (let index = 0; index < text.length; index++)
            hash = Math.imul(hash ^ text.charCodeAt(index), 16777619);
        return (hash >>> 0).toString(36);
    };
    const lengths = [0];
    for (let index = 1; index < points.length; index++)
        lengths.push(lengths[index - 1] + Math.hypot(points[index][0] - points[index - 1][0], points[index][1] - points[index - 1][1]));
    const at = (length) => {
        let index = 1;
        while (index < points.length - 1 && lengths[index] < length)
            index++;
        const span = lengths[index] - lengths[index - 1] || 1;
        const t = Math.max(0, Math.min(1, (length - lengths[index - 1]) / span));
        return [points[index - 1][0] + (points[index][0] - points[index - 1][0]) * t, points[index - 1][1] + (points[index][1] - points[index - 1][1]) * t];
    };
    const pieceOf = (from, to) => [at(from), ...points.filter((_, index) => lengths[index] > from + 1e-6 && lengths[index] < to - 1e-6), at(to)];
    const total = lengths[lengths.length - 1];
    if (total < 1e-6)
        return "";
    let [start, end] = Array.isArray(fade) ? fade : [fade, fade];
    if (start + end > total) {
        const share = total / (start + end);
        start *= share;
        end *= share;
    }
    const paint = `var(--anatomy-${tone})`;
    const ramp = (piece, rising) => {
        const d = pathOf(piece);
        const [x1, y1] = piece[0];
        const [x2, y2] = piece[piece.length - 1];
        const id = `iso-fade-${hashOf(`${d}|${tone}|${rising}`)}`;
        const marks = Array.from({ length: stops + 1 }, (_, index) => {
            const t = index / stops;
            const alpha = ease(rising ? t : 1 - t);
            return `<stop offset="${t.toFixed(3)}" stop-opacity="${alpha.toFixed(3)}" style="stop-color:${paint}"/>`;
        }).join("");
        const gradient = `<linearGradient id="${id}" gradientUnits="userSpaceOnUse" x1="${tenth(x1)}" y1="${tenth(y1)}" x2="${tenth(x2)}" y2="${tenth(y2)}">${marks}</linearGradient>`;
        return `<defs>${gradient}</defs><path${attrs({ class: className ? `iso-line ${className}` : "iso-line", "data-tone": tone, "data-free": rising ? "start" : "end", d, style: `stroke:url(#${id})` })}/>`;
    };
    const parts = [];
    if (start > 0.05)
        parts.push(ramp(pieceOf(0, start), true));
    if (total - end - start > 0.05)
        parts.push(lineSvg(pathOf(pieceOf(start, total - end)), { tone, className }));
    if (end > 0.05)
        parts.push(ramp(pieceOf(total - end, total), false));
    return `<g class="iso-fade">${parts.join("")}</g>`;
}
function wireSvg(from, to, { via = [], tone = "lo", dotted = false, dashed = false, flow = false, ends = "both", size = 1.2 } = {}) {
    const d = pathOf([from, ...via, to]);
    const end = ([x, y]) => `<circle class="iso-wire-end" cx="${tenth(x)}" cy="${tenth(y)}" r="${size}"/>`;
    return `<g class="iso-wire" data-tone="${tone}">${lineSvg(d, { tone, dotted, dashed, flow })}${ends === "both" || ends === "from" ? end(from) : ""}${ends === "both" || ends === "to" ? end(to) : ""}</g>`;
}
function faceTextSvg(transform, text, { size = 3, tone = "mid", anchor = "start" } = {}) {
    return `<text class="iso-face-text" data-tone="${tone}" transform="${transform}" font-size="${size}" text-anchor="${anchor}">${escape(text)}</text>`;
}
function groupSvg(children, extra = {}) {
    return `<g${attrs(extra)}>${Array.isArray(children) ? children.join("") : children}</g>`;
}
function figureSvg({ width, height, label, body, standalone = false, className }) {
    const style = standalone ? `<style>${ISO_CSS}</style>` : "";
    return `<svg${attrs({
        xmlns: "http://www.w3.org/2000/svg",
        class: `iso-svg${standalone ? " iso" : ""}${className ? ` ${className}` : ""}`,
        "data-theme": standalone || undefined,
        viewBox: `0 0 ${width} ${height}`,
        role: "img",
        "aria-label": label,
    })}>${style}${Array.isArray(body) ? body.join("") : body}</svg>`;
}
function swatchSvg(mark) {
    if (mark === "raised" || mark === "flat" || mark === "lit") {
        const body = mark === "flat" ? "M2 6.5 8 3.5 14 6.5 8 9.5Z" : "M2 4.5 8 1.5 14 4.5 14 6 8 9 2 6Z";
        const crease = mark === "flat" ? "" : `<path class="iso-swatch-crease" d="M2 4.5 8 7.5 14 4.5"/>`;
        return `<svg class="iso-swatch" data-mark="${mark}" viewBox="0 0 16 10" aria-hidden="true"><path d="${body}"/>${crease}</svg>`;
    }
    return `<svg class="iso-swatch" data-mark="${mark}" viewBox="0 0 16 10" aria-hidden="true"><path d="M1.5 5h13"/></svg>`;
}
function plateHtml({ fig, title, hint, readout = "", keys = [], caption, body }) {
    const corner = (name, text, extra = "") => `<span class="iso-plate-corner" data-corner="${name}"${extra}>${escape(text)}</span>`;
    const legend = keys.length || caption
        ? `<figcaption class="iso-legend">${keys.length ? `<span class="iso-keys">${keys.map((key) => `<span class="iso-key">${swatchSvg(key.mark)}${escape(key.label)}</span>`).join("")}</span>` : ""}${caption ? `<span class="iso-caption">${escape(caption)}</span>` : ""}</figcaption>`
        : "";
    return `<figure class="iso-figure"><div class="iso-plate"><div class="iso-plate-corners" aria-hidden="true">${corner("fig", fig)}${title ? corner("title", title) : ""}${hint ? corner("hint", hint) : ""}${corner("readout", readout, " data-readout")}</div>${body}</div>${legend}</figure>`;
}
function pageHtml({ title, theme = "dark", body, script = "", width = 650 }) {
    return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>${escape(title)}</title><style>${ISO_CSS}${PAGE_CSS.replace("__WIDTH__", String(width))}</style></head><body class="iso" data-theme="${theme}"><main class="iso-page">${body}</main>${script ? `<script type="module">${script.replace(/<\/(script)/gi, "<\\/$1")}<\/script>` : ""}<script>
(function(){
  const back = document.getElementById('btn-volver-calc');
  if (back && location.protocol === 'file:') back.href = 'index.html';
  const q = new URLSearchParams(location.search);
  if ([...q.keys()].length && window.BiorreactorApp) {
    const num = (ks) => { for (const k of ks) { const v = parseFloat(q.get(k)); if (Number.isFinite(v) && v > 0) return v; } return null; };
    const p = {}, Di = num(['Di','Dt']), Ht = num(['Ht','H']), De = num(['De','dt']), Dc = num(['Dc']), paso = num(['p']), Vw = num(['Vw','V']), N = num(['N']), Tin = num(['Tin','T_c_in']), Tout = num(['Tout','T_c_out']), caudal = num(['caudal','Q_agua']), T0 = num(['T0','T_op']);
    if (Di) p.Di = Di < 5 ? Di * 100 : Di;
    if (Ht) p.Ht = Ht < 10 ? Ht * 100 : Ht;
    if (De) p.De = De < 0.2 ? De * 1000 : De;
    if (Dc) p.Dc = Dc < 5 ? Dc * 100 : Dc;
    if (paso) p.p = paso < 1 ? paso * 100 : paso;
    if (Vw) p.Vw = Vw < 5 ? Vw * 1000 : Vw;
    if (N) p.N = N;
    if (Tin !== null) p.Tin = Tin;
    if (Tout !== null) p.Tout = Tout;
    if (caudal) p.caudal = caudal;
    if (T0) { p.T0 = T0; p.Tobj = Math.max(4, T0 - 10); }
    try { window.BiorreactorApp.setParametros(p); } catch (e) {}
  }
})();
</script></body></html>`;
}
const PAGE_CSS = `
*{box-sizing:border-box}
html{color-scheme:dark}
body[data-theme="light"]{color-scheme:light}
body{margin:0;background:var(--anatomy-page);color:var(--anatomy-text);font:14px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;-webkit-font-smoothing:antialiased}
.iso-page{width:min(__WIDTH__px,100% - 32px);margin:64px auto 96px}
`;
const ISO_CSS = `
.iso{
  --anatomy-page:#0b0b0b;--anatomy-card:#131313;--anatomy-text:#a3a3a3;--anatomy-strong:#e6e6e6;--anatomy-muted:#6e6e6e;
  --anatomy-paper:#131313;--anatomy-face:#1b1b1b;--anatomy-top:#1c1c1c;
  --anatomy-shade-0:#0f0f0f;--anatomy-shade-1:#121212;--anatomy-shade-2:#151515;--anatomy-shade-3:#181818;
  --anatomy-lit-top:#262626;--anatomy-lit-shade:#1d1d1d;
  --anatomy-hi:#d4d4d4;--anatomy-mid:#6a6a6a;--anatomy-lo:#3e3e3e;--anatomy-faint:#2a2a2a;--anatomy-lit:#ffffff;
  --anatomy-halo:#050505;--anatomy-halo-opacity:.9;--anatomy-dot:#e6e6e6;
  --anatomy-red:#ff6b6b;--anatomy-green:#6be0a0;--anatomy-blue:#6b9bff;
  --anatomy-weight:.6px;--anatomy-ease:cubic-bezier(.32,.72,0,1);
  --anatomy-mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  --anatomy-lift:inset 0 1px 0 0 #ffffff0a,inset 0 0 0 1px #ffffff06,0 0 0 1px #00000026,0 1px 1px -.5px #00000014,0 3px 3px -1.5px #00000014,0 6px 6px -3px #00000014,0 12px 12px -6px #00000014;
}
.iso[data-theme="light"]{
  --anatomy-page:#f2f2ef;--anatomy-card:#fbfbfa;--anatomy-text:#57575a;--anatomy-strong:#18181a;--anatomy-muted:#8c8c90;
  --anatomy-paper:#fbfbfa;--anatomy-face:#ffffff;--anatomy-top:#ffffff;
  --anatomy-shade-0:#dcdcd9;--anatomy-shade-1:#e5e5e2;--anatomy-shade-2:#ededea;--anatomy-shade-3:#f4f4f1;
  --anatomy-lit-top:#ffffff;--anatomy-lit-shade:#eeeeeb;
  --anatomy-hi:#26262a;--anatomy-mid:#8a8a8f;--anatomy-lo:#b4b4b8;--anatomy-faint:#d4d4d6;--anatomy-lit:#000000;
  --anatomy-halo:#8a8a80;--anatomy-halo-opacity:.32;--anatomy-dot:#3a3a3e;
  --anatomy-red:#e5484d;--anatomy-green:#30a46c;--anatomy-blue:#3e63dd;
  --anatomy-lift:inset 0 1px 0 0 #ffffff,0 0 0 1px #0000000f,0 1px 1px -.5px #0000000a,0 3px 3px -1.5px #0000000a,0 6px 6px -3px #0000000a,0 12px 12px -6px #0000000a;
}
.iso-svg{display:block;width:100%;height:auto;overflow:visible}
.iso-svg path,.iso-svg circle{vector-effect:non-scaling-stroke}
.iso-fill{fill:var(--anatomy-paper);stroke:none}
.iso-solid[data-lit]>.iso-fill{fill:var(--anatomy-face)}
.iso-shade,.iso-top{stroke:none}
.iso-shade[data-shade="0"]{fill:var(--anatomy-shade-0)}
.iso-shade[data-shade="1"]{fill:var(--anatomy-shade-1)}
.iso-shade[data-shade="2"]{fill:var(--anatomy-shade-2)}
.iso-shade[data-shade="3"]{fill:var(--anatomy-shade-3)}
.iso-top{fill:var(--anatomy-top)}
.iso-solid[data-lit]>.iso-top{fill:var(--anatomy-lit-top)}
.iso-solid[data-lit]>.iso-shade{fill:var(--anatomy-lit-shade)}
.iso-solid[data-seamless]>.iso-shade,.iso-solid[data-seamless]>.iso-top{stroke-width:.4;stroke-linejoin:round}
.iso-solid[data-seamless]>.iso-shade[data-shade="0"]{stroke:var(--anatomy-shade-0)}
.iso-solid[data-seamless]>.iso-shade[data-shade="1"]{stroke:var(--anatomy-shade-1)}
.iso-solid[data-seamless]>.iso-shade[data-shade="2"]{stroke:var(--anatomy-shade-2)}
.iso-solid[data-seamless]>.iso-shade[data-shade="3"]{stroke:var(--anatomy-shade-3)}
.iso-solid[data-seamless]>.iso-top{stroke:var(--anatomy-top)}
.iso-solid[data-seamless][data-lit]>.iso-top{stroke:var(--anatomy-lit-top)}
.iso-solid[data-seamless][data-lit]>.iso-shade{stroke:var(--anatomy-lit-shade)}
.iso-solid[data-lit]>.iso-edge{stroke:var(--anatomy-lit)}
.iso-halo{fill:var(--anatomy-halo);stroke:none;filter:blur(9px);opacity:var(--anatomy-halo-opacity);transform:translateY(3px)}
.iso-line{fill:none;stroke:var(--anatomy-lo);stroke-width:var(--anatomy-weight);stroke-linecap:round;stroke-linejoin:round}
.iso-line[data-tone="hi"]{stroke:var(--anatomy-hi)}
.iso-line[data-tone="mid"]{stroke:var(--anatomy-mid)}
.iso-line[data-tone="lo"]{stroke:var(--anatomy-lo)}
.iso-line[data-tone="faint"]{stroke:var(--anatomy-faint)}
.iso-line[data-tone="lit"]{stroke:var(--anatomy-lit)}
.iso-line[data-tone="red"]{stroke:var(--anatomy-red)}
.iso-line[data-tone="green"]{stroke:var(--anatomy-green)}
.iso-line[data-tone="blue"]{stroke:var(--anatomy-blue)}
.iso-line[data-dotted]{stroke-dasharray:.01 3.5;stroke-width:1}
.iso-line[data-dashed]{stroke-dasharray:3 3}
.iso-line[data-flow]{stroke-dasharray:2 6}
.iso-wire-end{fill:var(--anatomy-paper);stroke:var(--anatomy-lo);stroke-width:var(--anatomy-weight)}
.iso-wire[data-tone="hi"] .iso-wire-end{stroke:var(--anatomy-hi)}
.iso-wire[data-tone="mid"] .iso-wire-end{stroke:var(--anatomy-mid)}
.iso-wire[data-tone="faint"] .iso-wire-end{stroke:var(--anatomy-faint)}
.iso-wire[data-tone="lit"] .iso-wire-end{stroke:var(--anatomy-lit)}
.iso-dots circle{fill:var(--anatomy-dot);stroke:none;opacity:.32}
.iso-dots[data-tone="hi"] circle{opacity:.85}
.iso-dots[data-tone="mid"] circle{opacity:.45}
.iso-dots[data-tone="lo"] circle{opacity:.22}
.iso-dots[data-tone="faint"] circle{opacity:.12}
.iso-face-text{fill:var(--anatomy-mid);font-family:var(--anatomy-mono);letter-spacing:.04em}
.iso-face-text[data-tone="hi"]{fill:var(--anatomy-hi)}
.iso-face-text[data-tone="lo"]{fill:var(--anatomy-lo)}
.iso-face-text[data-tone="faint"]{fill:var(--anatomy-faint)}
.tb-body{fill:var(--anatomy-shade-2);stroke:none}
.tb-shine{fill:var(--anatomy-top);stroke:none}
.tb-shade{fill:var(--anatomy-shade-0);stroke:none}
.tb-ring{fill:none;stroke:var(--anatomy-lo);stroke-width:.45;stroke-linecap:butt;stroke-linejoin:round}
.tb-ring[data-tone="mid"]{stroke:var(--anatomy-mid)}
.tb-ring[data-tone="faint"]{stroke:var(--anatomy-faint)}
.tb-ring[data-tone="hi"]{stroke:var(--anatomy-hi)}
@keyframes iso-flow{to{stroke-dashoffset:-32}}
@keyframes iso-pulse{0%,100%{opacity:.18}12%{opacity:1}40%{opacity:.18}}
.iso-line[data-flow]{animation:iso-flow 1600ms linear infinite}
.iso-dots circle[data-pulse]{animation:iso-pulse 1600ms steps(1,end) infinite;animation-delay:var(--pulse-at,0ms)}
.iso-figure{margin:32px 0 36px}
.iso-plate{position:relative;padding:52px 24px;border-radius:0;background:var(--anatomy-card);box-shadow:0 0 0 1px var(--anatomy-mid);--anatomy-paper:var(--anatomy-card)}
.iso-plate-corners{position:absolute;inset:18px 20px;pointer-events:none}
.iso-plate-corner{position:absolute;color:var(--anatomy-muted);font-family:var(--anatomy-mono);font-size:10.5px;line-height:1;letter-spacing:.04em;white-space:nowrap}
.iso-plate-corner[data-corner="fig"]{top:0;left:0;color:var(--anatomy-text)}
.iso-plate-corner[data-corner="title"]{top:0;right:0;text-transform:uppercase;letter-spacing:.1em}
.iso-plate-corner[data-corner="hint"]{bottom:0;left:0;text-transform:uppercase;letter-spacing:.1em}
.iso-plate-corner[data-corner="readout"]{right:0;bottom:0;color:var(--anatomy-text);font-variant-numeric:tabular-nums}
.iso-plate .iso-svg{margin:0 auto}
.iso-legend{display:grid;gap:8px;margin-top:12px;color:var(--anatomy-muted);font-size:13px;line-height:1.55}
.iso-keys{display:flex;flex-wrap:wrap;gap:6px 18px}
.iso-key{display:inline-flex;align-items:center;gap:7px}
.iso-swatch{width:16px;height:10px;flex:none;overflow:visible}
.iso-swatch path{fill:none;stroke:var(--anatomy-hi);stroke-width:1;vector-effect:non-scaling-stroke;stroke-linejoin:round;stroke-linecap:round}
.iso-swatch[data-mark="raised"] path,.iso-swatch[data-mark="flat"] path{stroke:var(--anatomy-mid)}
.iso-swatch[data-mark="lit"] path{stroke:var(--anatomy-lit)}
.iso-swatch[data-mark="dim"] path{stroke:var(--anatomy-lo)}
.iso-swatch[data-mark="dotted"] path{stroke-dasharray:.01 3;stroke-width:1.4}
.iso-swatch[data-mark="dashed"] path{stroke-dasharray:3 2.5}
.iso-swatch[data-mark="red"] path{stroke:var(--anatomy-red)}
.iso-swatch[data-mark="green"] path{stroke:var(--anatomy-green)}
.iso-swatch[data-mark="blue"] path{stroke:var(--anatomy-blue)}
.iso-swatch-crease{stroke:var(--anatomy-lo)!important}
@media (max-width:560px){.iso-plate{padding:48px 12px 70px}.iso-plate-corner[data-corner="hint"]{bottom:18px}}
@media (prefers-reduced-motion:reduce){.iso-line[data-flow],.iso-dots circle[data-pulse]{animation:none}}
`;

// Cálculos de ingeniería del serpentín helicoidal de un biorreactor.
// Módulo puro (sin DOM). Unidades de entrada: las de la interfaz (ver CAMPOS).
// Todos los resultados son aproximados y se basan en los supuestos de SUPUESTOS.

const SUPUESTOS = {
  rhoMedio: 1000, // kg/m3, medio de cultivo tratado como agua
  cpMedio: 4186, // J/(kg·K)
  rhoRef: 1.0, // kg/L, refrigerante tratado como agua
  cpRef: 4180, // J/(kg·K)
  tolerancia: 0.1, // diferencia relativa aceptada entre los dos balances
};

// Valores ILUSTRATIVOS para poder dibujar el modelo. No son resultados experimentales.
const VALORES_ILUSTRATIVOS = {
  Vw: 50, // L
  Di: 40, // cm
  Ht: 60, // cm
  De: 12.7, // mm
  Dc: 22, // cm
  N: 8, // vueltas
  p: 3, // cm
  T0: 37, // °C
  Tobj: 25, // °C
  Tin: 8, // °C
  Tout: 12, // °C
  caudal: 11, // L/min
  U: 500, // W/(m2·K)
};

const CAMPOS = {
  Vw: { etiqueta: "Volumen de trabajo", unidad: "L", min: 1, max: 2000, paso: 1, grupo: "geo" },
  Di: { etiqueta: "Diámetro interno del tanque", unidad: "cm", min: 15, max: 150, paso: 0.5, grupo: "geo" },
  Ht: { etiqueta: "Altura del tanque", unidad: "cm", min: 25, max: 300, paso: 0.5, grupo: "geo" },
  De: { etiqueta: "Diámetro exterior de la tubería", unidad: "mm", min: 3, max: 50, paso: 0.1, grupo: "geo" },
  Dc: { etiqueta: "Diámetro de la hélice", unidad: "cm", min: 5, max: 140, paso: 0.5, grupo: "geo" },
  L: { etiqueta: "Longitud del serpentín", unidad: "m", min: 0.1, max: 200, paso: 0.05, grupo: "geo" },
  N: { etiqueta: "Número de vueltas", unidad: "vueltas", min: 1, max: 40, paso: 0.25, grupo: "geo" },
  p: { etiqueta: "Paso entre vueltas", unidad: "cm", min: 0.5, max: 30, paso: 0.1, grupo: "geo" },
  T0: { etiqueta: "Temperatura inicial del cultivo", unidad: "°C", min: -5, max: 80, paso: 0.5, grupo: "term" },
  Tobj: { etiqueta: "Temperatura objetivo del cultivo", unidad: "°C", min: -5, max: 80, paso: 0.5, grupo: "term" },
  Tin: { etiqueta: "Temperatura de entrada del refrigerante", unidad: "°C", min: -20, max: 60, paso: 0.5, grupo: "term" },
  Tout: { etiqueta: "Temperatura de salida del refrigerante", unidad: "°C", min: -20, max: 80, paso: 0.5, grupo: "term" },
  caudal: { etiqueta: "Caudal del refrigerante", unidad: "L/min", min: 0.1, max: 500, paso: 0.1, grupo: "term" },
  U: { etiqueta: "Coeficiente global de transferencia", unidad: "W/(m²·K)", min: 10, max: 5000, paso: 10, grupo: "term" },
};

// L = N · √[(π·Dc)² + p²]   (Dc y p en la misma unidad que L)
function longitudHelice(N, Dc, p) {
  return N * Math.sqrt((Math.PI * Dc) ** 2 + p ** 2);
}

// Inversa: N = L / √[(π·Dc)² + p²]
function vueltasDesdeLongitud(L, Dc, p) {
  return L / Math.sqrt((Math.PI * Dc) ** 2 + p ** 2);
}

// A = π · De · L   (superficie exterior de la tubería, sin considerar extremos ni soportes)
function areaSerpentin(De, L) {
  return Math.PI * De * L;
}

// ΔTlm = (ΔT1 − ΔT2) / ln(ΔT1/ΔT2); si ΔT1 = ΔT2 el límite es ΔT1.
function dtMediaLog(dT1, dT2) {
  if (!(dT1 > 0) || !(dT2 > 0)) return NaN;
  if (Math.abs(dT1 - dT2) < 1e-9 * Math.max(dT1, dT2)) return dT1;
  return (dT1 - dT2) / Math.log(dT1 / dT2);
}

// Completa L a partir de N (o N a partir de L) según el campo editado.
function sincronizar(par, editado) {
  const q = { ...par };
  const Dc = q.Dc / 100;
  const p = q.p / 100;
  if (editado === "L") q.N = vueltasDesdeLongitud(q.L, Dc, p);
  else q.L = longitudHelice(q.N, Dc, p);
  return q;
}

const esNum = (v) => typeof v === "number" && Number.isFinite(v);

// Geometría derivada. Las proporciones del impulsor y los deflectores son ILUSTRATIVAS
// (turbina Rushton Di/3, deflectores Di/10); no son datos del usuario.
function geometria(par) {
  const Di = par.Di / 100;
  const Ht = par.Ht / 100;
  const De = par.De / 1000;
  const Dc = par.Dc / 100;
  const p = par.p / 100;
  const Vw = par.Vw / 1000;
  const Hl = (4 * Vw) / (Math.PI * Di * Di);
  const Dimp = Di / 3;
  const Wb = Di / 10;
  const Rc = Dc / 2;
  const Rr = Rc + 1.8 * De; // radio del tramo de entrada (por fuera de la hélice)
  const z0 = Math.max(0.08 * Di, De); // base de la hélice sobre el fondo
  const altura = par.N * p;
  return { Di, Ht, De, Dc, p, Vw, Hl, Dimp, Wb, Rc, Rr, z0, altura, zTop: z0 + altura, holgura: 0.01 };
}

function validar(par) {
  const errores = [];
  const avisos = [];
  const err = (campo, msg) => errores.push({ campo, msg });
  const avi = (campo, msg) => avisos.push({ campo, msg });

  for (const [id, c] of Object.entries(CAMPOS)) {
    const v = par[id];
    if (!esNum(v)) err(id, `${c.etiqueta}: escribe un número.`);
    else if (v < c.min || v > c.max) err(id, `${c.etiqueta}: usa un valor entre ${c.min} y ${c.max} ${c.unidad}.`);
  }
  if (errores.length) return { errores, avisos, geo: null };

  const g = geometria(par);
  const mm = (m) => (m * 1000).toFixed(1);
  const cm = (m) => (m * 100).toFixed(1);

  if (g.Hl > g.Ht) err("Vw", `El volumen no cabe: llenaría ${cm(g.Hl)} cm y el tanque mide ${cm(g.Ht)} cm.`);
  else if (g.Hl / g.Ht > 0.8) avi("Vw", "El nivel supera el 80 % de la altura: queda poco espacio de cabeza.");
  if (g.Dc <= g.De) err("Dc", "El diámetro de la hélice debe ser mayor que el de la tubería.");
  if (g.p < g.De) err("p", `El paso (${cm(g.p)} cm) es menor que el diámetro exterior (${mm(g.De)} mm): las vueltas se solapan.`);
  else if (g.p < 1.5 * g.De) avi("p", "Paso menor que 1,5 veces el diámetro exterior: poco espacio de circulación entre espiras.");

  const dcMin = g.Dimp + g.De + 2 * g.holgura;
  const dcMax = 2 * (g.Di / 2 - g.Wb - g.holgura - 2.3 * g.De);
  if (g.Dc < dcMin) err("Dc", `La hélice toca el impulsor (Di/3 = ${cm(g.Dimp)} cm). Usa un diámetro de hélice mayor que ${cm(dcMin)} cm.`);
  if (g.Dc > dcMax) err("Dc", `La hélice (con su tramo de entrada) choca con los deflectores. Usa un diámetro de hélice menor que ${cm(dcMax)} cm.`);
  if (g.zTop + g.De / 2 > g.Hl) err("N", `El serpentín mide ${cm(g.zTop + g.De / 2)} cm de alto y el nivel es ${cm(g.Hl)} cm: queda fuera del medio.`);

  if (par.T0 <= par.Tobj) err("T0", "La temperatura inicial debe ser mayor que la objetivo para que haya enfriamiento.");
  if (par.Tin >= par.Tout) err("Tout", "El refrigerante se calienta al pasar: la salida debe ser mayor que la entrada.");
  if (par.Tout >= par.Tobj) err("Tout", "La salida del refrigerante debe ser menor que la temperatura objetivo del cultivo.");
  if (par.U < 100 || par.U > 3000) avi("U", "Valor poco habitual para agua-agua con serpentín (orientativo). Conviene verificarlo con una correlación.");
  return { errores, avisos, geo: g };
}

// Calcula todo a la temperatura de cultivo Tc (por defecto la inicial).
function calcular(par, Tc) {
  const v = validar(par);
  if (v.errores.length) return { ok: false, errores: v.errores, avisos: v.avisos };
  const g = v.geo;
  const S = SUPUESTOS;
  const T = esNum(Tc) ? Math.min(par.T0, Math.max(par.Tobj, Tc)) : par.T0;

  const masa = S.rhoMedio * g.Vw; // kg
  const energia = masa * S.cpMedio * (par.T0 - par.Tobj); // J

  const L = longitudHelice(par.N, g.Dc, g.p);
  const A = areaSerpentin(g.De, L);
  const gasto = (par.caudal * S.rhoRef) / 60; // kg/s
  const qRef = gasto * S.cpRef * (par.Tout - par.Tin); // W

  const en = (Tcult) => {
    const dT1 = Tcult - par.Tin;
    const dT2 = Tcult - par.Tout;
    const dtlm = dtMediaLog(dT1, dT2);
    const q = par.U * A * dtlm;
    return { Tcult, dT1, dT2, dtlm, q };
  };
  const ini = en(par.T0);
  const obj = en(par.Tobj);
  const act = en(T);
  const dif = Math.abs(qRef - act.q) / act.q;
  const estado = dif <= S.tolerancia ? "coherente" : dif <= 2.5 * S.tolerancia ? "revisar" : "inconsistente";
  const areaNecesaria = qRef / (par.U * act.dtlm);

  return {
    ok: true,
    avisos: v.avisos,
    geo: g,
    T,
    masa,
    energia,
    L,
    A,
    gasto,
    qRef,
    ini,
    obj,
    act,
    dif,
    estado,
    areaNecesaria,
    llenado: g.Hl / g.Ht,
    volTanque: (Math.PI / 4) * g.Di * g.Di * g.Ht * 1000, // L
  };
}

// Textos de los componentes. Los valores numéricos de "params" se leen de los parámetros actuales.
const COMPONENTES = {
  shell: {
    nombre: "Cuerpo del tanque",
    funcion: "Contener el medio de cultivo y resistir la presión de operación y de esterilización.",
    descripcion: "Cilindro vertical de pared lisa y fondo plano, con soldaduras circunferenciales y acabado interior pulido para facilitar la limpieza. En la ilustración se retira una cuña frontal para ver el interior; la vista Transparente y la Cerrada se eligen arriba.",
    material: "Acero inoxidable 316L con pulido interior (sugerido; depende del proceso y de la norma aplicable).",
    importancia: "Su diámetro y su altura fijan el volumen disponible y el espacio donde cabe el serpentín.",
    params: ["Di", "Ht", "Vw"],
  },
  lid: {
    nombre: "Tapa superior",
    funcion: "Cerrar el tanque y servir de apoyo al motor, a los sensores, a los puertos y a las conexiones del serpentín.",
    descripcion: "Placa circular embridada con junta y tornillos perimetrales. Tiene bocas para inoculación, muestreo, filtro de gases de salida, sensores y los dos pasos del refrigerante.",
    material: "Acero inoxidable 316L; juntas de elastómero apto para esterilización (EPDM o silicona, según el proceso).",
    importancia: "Mantiene la barrera aséptica. Los soportes del serpentín cuelgan de ella.",
    params: ["Di"],
  },
  motor: {
    nombre: "Motor y reductor",
    funcion: "Hacer girar el eje de agitación a la velocidad de proceso.",
    descripcion: "Motor eléctrico con reductor montado sobre la tapa. Transmite el giro al eje a través de un sello mecánico. La velocidad de giro no es un parámetro de este modelo y la animación es solo ilustrativa.",
    material: "Carcasa de aluminio o acero con protección contra salpicaduras; piezas en contacto con el proceso en 316L.",
    importancia: "La agitación acerca el calor del medio a la pared del serpentín. Con poca mezcla, el coeficiente U real baja.",
    params: [],
  },
  shaft: {
    nombre: "Eje de agitación",
    funcion: "Transmitir el par del motor a los impulsores.",
    descripcion: "Barra cilíndrica que atraviesa la tapa por el sello mecánico y termina en un apoyo cerca del fondo. Los impulsores van fijos en él.",
    material: "Acero inoxidable 316L, superficie pulida.",
    importancia: "Su diámetro y su longitud limitan la velocidad y el par que se puede aplicar sin vibración.",
    params: ["Ht"],
  },
  impeller: {
    nombre: "Impulsores (turbina Rushton)",
    funcion: "Mezclar el medio y generar un flujo radial que reparte el calor y los nutrientes.",
    descripcion: "Disco con seis palas planas verticales. Producen flujo radial hacia la pared, que se divide en dos lazos de recirculación, uno por encima y otro por debajo. En la ilustración el diámetro del impulsor es Di/3, una proporción de referencia y no un dato del usuario.",
    material: "Acero inoxidable 316L.",
    importancia: "El serpentín debe quedar fuera del barrido de las palas, por eso la hélice no puede ser más estrecha que el impulsor.",
    params: ["Di"],
  },
  baffle: {
    nombre: "Deflectores",
    funcion: "Evitar el vórtice central y convertir el giro del líquido en mezcla vertical.",
    descripcion: "Cuatro placas verticales pegadas a la pared a 90° entre sí. Su ancho en la ilustración es Di/10, una proporción de referencia. Ocupan la zona exterior, así que limitan el diámetro máximo de la hélice.",
    material: "Acero inoxidable 316L.",
    importancia: "Sin deflectores el líquido gira como un bloque y el intercambio de calor con el serpentín empeora.",
    params: ["Di"],
  },
  medium: {
    nombre: "Medio de cultivo",
    funcion: "Alojar a las células y aportar los nutrientes. Es la fuente del calor que hay que retirar.",
    descripcion: "En los cálculos se trata como agua (densidad 1000 kg/m³, Cp 4186 J/(kg·K)) y como mezcla uniforme a una sola temperatura. El color del medio sigue la temperatura mostrada, y la zona azulada alrededor del serpentín es solo una representación del enfriamiento local.",
    material: "No aplica. Líquido acuoso con sales, fuente de carbono y nitrógeno.",
    importancia: "Su masa y su salto de temperatura fijan la energía que el serpentín debe retirar (Q = m·Cp·ΔT).",
    params: ["Vw", "T0", "Tobj"],
  },
  tprobe: {
    nombre: "Sensor de temperatura",
    funcion: "Medir la temperatura del cultivo para el lazo que regula el caudal del refrigerante.",
    descripcion: "Sonda Pt100 dentro de una vaina que se introduce por un prensaestopas de la tapa. Queda lejos del serpentín para no medir el agua fría de la tubería.",
    material: "Elemento Pt100 en vaina de acero inoxidable 316L.",
    importancia: "El control de temperatura depende de esta lectura.",
    params: ["Tobj"],
  },
  phprobe: {
    nombre: "Sensor de pH",
    funcion: "Medir el pH del medio.",
    descripcion: "Electrodo de vidrio combinado con referencia, esterilizable, montado en un puerto de la tapa. Se calibra antes de esterilizar.",
    material: "Cuerpo de vidrio o polímero apto para esterilización; conector fuera del tanque.",
    importancia: "El pH condiciona el crecimiento y suele controlarse con adición de ácido o base.",
    params: [],
  },
  inlet: {
    nombre: "Entrada del medio",
    funcion: "Alimentar medio fresco o reactivos al tanque.",
    descripcion: "Boquilla lateral con conexión sanitaria. En la ilustración está en la zona alta de la pared.",
    material: "Acero inoxidable 316L con unión sanitaria y junta de elastómero.",
    importancia: "Permite operar en lote alimentado o en continuo.",
    params: ["Vw"],
  },
  outlet: {
    nombre: "Salida del medio",
    funcion: "Retirar medio para cosecha, muestreo o recirculación.",
    descripcion: "Boquilla lateral cerca del fondo con conexión sanitaria, en el lado opuesto a la entrada.",
    material: "Acero inoxidable 316L con unión sanitaria.",
    importancia: "Debe quedar por debajo del nivel de operación para poder vaciar el tanque.",
    params: ["Vw"],
  },
  drain: {
    nombre: "Drenaje de fondo",
    funcion: "Vaciar el tanque por completo al terminar el lote o durante la limpieza.",
    descripcion: "Tubo central bajo el fondo con una válvula de bola sanitaria y una conexión de salida. Va montado a ras del fondo para evitar bolsas de líquido.",
    material: "Válvula y tubería de acero inoxidable 316L, asientos de PTFE.",
    importancia: "Un drenaje con zonas muertas deja residuos que contaminan el siguiente lote.",
    params: [],
  },
  coil: {
    nombre: "Serpentín helicoidal",
    funcion: "Transferir calor entre el medio de cultivo y el refrigerante que circula por dentro de la tubería.",
    descripcion: "Tubería enrollada en hélice alrededor del eje, entre el impulsor y los deflectores. El refrigerante entra por el tubo exterior, sube por la hélice y sale por la tapa. Longitud aproximada L = N·√[(π·Dc)² + p²] y superficie A = π·De·L; son aproximaciones y no un modelo térmico completo.",
    material: "Acero inoxidable 316L, sujeto a los requisitos del proceso.",
    importancia: "Su superficie y el coeficiente U determinan cuánto calor se puede retirar con una diferencia de temperatura dada.",
    params: ["De", "Dc", "N", "p", "L"],
  },
  cin: {
    nombre: "Entrada del refrigerante",
    funcion: "Recibir el refrigerante frío desde la red de agua helada o el equipo de enfriamiento.",
    descripcion: "El tubo baja por la tapa, recorre el fondo por fuera de la hélice y entra por la parte inferior del serpentín. La conexión exterior es una unión sanitaria o de manguera.",
    material: "Acero inoxidable 316L; manguera flexible apta para la temperatura de trabajo.",
    importancia: "La temperatura y el caudal de entrada fijan cuánto calor puede absorber el refrigerante.",
    params: ["Tin", "caudal"],
  },
  cout: {
    nombre: "Salida del refrigerante",
    funcion: "Devolver el refrigerante ya calentado al equipo de enfriamiento.",
    descripcion: "El tubo sale por la parte alta de la hélice y atraviesa la tapa. Su temperatura es una de las entradas del balance del refrigerante.",
    material: "Acero inoxidable 316L; manguera flexible.",
    importancia: "Si la temperatura de salida no es coherente con el caudal y el área, los dos balances de la pestaña Cálculos no coinciden.",
    params: ["Tout", "caudal"],
  },
  support: {
    nombre: "Soportes del serpentín",
    funcion: "Sujetar la hélice para que no vibre ni roce los impulsores.",
    descripcion: "Dos varillas verticales colgadas de la tapa, con abrazaderas que sujetan algunas espiras. El tubo de entrada hace de tercer apoyo.",
    material: "Acero inoxidable 316L.",
    importancia: "Una hélice sin sujetar se mueve con el flujo y puede fatigarse en la soldadura de los pasos de tapa.",
    params: ["N", "p"],
  },
};

let ORDEN_COMPONENTES = ["coil", "cin", "cout", "support", "impeller", "shaft", "baffle", "medium", "shell", "lid", "motor", "tprobe", "phprobe", "inlet", "outlet", "drain"];

const VARIANT = {
  id: "blueprint", scene: "sceneC.js", theme: "light",
  titulo: "Plano en sección del biorreactor", sub: "Alzado con cotas que siguen los parámetros", fig: "Lámina 4", corner: "Sección A-A",
  defaults: { Vw: 200, Di: 60, Ht: 100, De: 19, Dc: 32, N: 10, p: 5, T0: 36, Tobj: 24, Tin: 6, Tout: 11, caudal: 26, U: 550 },
  vistas: { cut: "Con cotas", glass: "Sin cotas", closed: "Solo serpentín" },
  css: `body{--bp-paper:#0c2a4d;--bp-line:#d8ecff;--bp-dim:#7fd0ff;--bp-grid:#16406f;--bp-fill:#12365f;--acc:#ffd166}body[data-theme="light"]{--bp-paper:#fbfbf8;--bp-line:#1d2b3a;--bp-dim:#0b6bcb;--bp-grid:#e6e9ee;--bp-fill:#e9eef5;--acc:#d9480f}.iso-plate.stagecard{background:var(--bp-paper);border-color:var(--bp-dim)}.iso-plate-corner{color:var(--bp-dim)!important}`,
};
Object.assign(COMPONENTES.shell, { descripcion: "En el plano se ve la pared cortada con rayado diagonal, que es la forma habitual de marcar material cortado. Las cotas de diámetro interno, altura y nivel siguen los valores del panel." });

// ===== Plano en sección (alzado, 2D) =====
const SCENE = (() => {
  const W = 620; let H = 760;
  const D2R = Math.PI / 180;
  const f1 = (n) => (Math.round(n * 10) / 10).toString();
  const pth = (pts, close) => pts.map((p, i) => (i ? "L" : "M") + f1(p[0]) + " " + f1(p[1])).join("") + (close ? "Z" : "");
  const MIX = (a, b, t) => a.map((v, i) => Math.round(v + (b[i] - v) * t));
  const STOPS = [[0, [37, 99, 235]], [0.3, [56, 189, 248]], [0.55, [253, 230, 138]], [0.8, [251, 146, 60]], [1, [239, 68, 68]]];
  function tempColor(u) {
    const x = Math.min(1, Math.max(0, u));
    for (let i = 1; i < STOPS.length; i++) if (x <= STOPS[i][0]) { const [x0, c0] = STOPS[i - 1], [x1, c1] = STOPS[i]; const c = MIX(c0, c1, (x - x0) / (x1 - x0)); return `rgb(${c})`; }
    return "rgb(239,68,68)";
  }
  const part = (id, inner) => `<g data-part="${id}">${inner}</g>`;
  const hit = (d, w) => `<path class="bp-hit" d="${d}" style="stroke-width:${w}"/>`;
  const rect = (x, y, w, h, cls = "bp-o") => `<rect class="${cls}" x="${f1(x)}" y="${f1(y)}" width="${f1(w)}" height="${f1(h)}"/>`;

  function model(par, compact) {
    const padX = compact ? 40 : 128, padT = 70, padB = 52;
    H = compact ? 760 : 640;
    const Hmm = par.Ht + 0.52 * par.Di + 0.26 * par.Di; // tanque + motor + patas
    const s = Math.min((W - 2 * padX) / (par.Di * 1.28), (H - padT - padB) / Hmm);
    const cx = W / 2, yG = H - padB;
    const g = { s, cx, yG };
    g.Ri = (par.Di / 2) * s; g.t = Math.max(3.2, 0.03 * par.Di * s); g.Ro = g.Ri + g.t;
    g.leg = 0.2 * par.Di * s; g.tb = Math.max(3.4, 0.03 * par.Di * s);
    g.yB = yG - g.leg; // base interior del fondo
    g.Hh = par.Ht * s; g.yT = g.yB - g.Hh; g.Hl = Math.min(g.Hh * 0.97, ((4000 * par.Vw) / (Math.PI * par.Di * par.Di)) * s);
    g.yL = g.yB - g.Hl;
    g.de = Math.max(3, (par.De / 10) * s); g.Rc = (par.Dc / 2) * s; g.pv = Math.max(g.de * 1.05, par.p * s);
    g.z0 = Math.max(0.08 * par.Di * s, g.de); g.N = par.N; g.Rr = g.Rc + 1.8 * g.de;
    g.Dimp = (par.Di / 3) * s; g.Rimp = g.Dimp / 2; g.Wb = (par.Di / 10) * s;
    g.lidT = Math.max(4, 0.03 * par.Di * s);
    return g;
  }
  const Y = (g, z) => g.yB - z; // altura sobre el fondo -> y de pantalla

  function build(par, opt) {
    const g = model(par, opt.compact);
    const mode = opt.view, dims = mode === "cut", only = mode === "closed";
    const { cx } = g;
    const Tin = par.Tin, Tout = par.Tout, T0 = par.T0;
    const cScale = (T) => tempColor((T - Tin) / Math.max(1, T0 - Tin));
    const out = {}, anchors = {};
    const xl = cx - g.Ri, xr = cx + g.Ri;

    // fondo cuadriculado y suelo
    let base = `<g class="bp-grid">`;
    for (let x = 0; x <= W; x += 20) base += `<path d="M${x} 0V${H}"/>`;
    for (let y = 0; y <= H; y += 20) base += `<path d="M0 ${y}H${W}"/>`;
    base += `</g><path class="bp-o" d="M${cx - g.Ro * 1.7} ${g.yG}H${cx + g.Ro * 1.7}"/>`;
    for (let x = cx - g.Ro * 1.7; x < cx + g.Ro * 1.7; x += 9) base += `<path class="bp-hatch" d="M${f1(x)} ${g.yG}l-6 7"/>`;
    // patas y placa del fondo
    let body = "";
    for (const sx of [-1, 1]) body += rect(cx + sx * g.Ro * 0.78 - 5, g.yB + g.tb, 10, g.leg - g.tb, "bp-o bp-f") + rect(cx + sx * g.Ro * 0.78 - 9, g.yG - 3, 18, 3, "bp-o bp-f");
    body += rect(cx - g.Ro - 3, g.yB, 2 * g.Ro + 6, g.tb, "bp-o bp-cut");
    if (!only) base += part("shell", body); else base += `<g class="bp-dim-only">${body}</g>`;
    out.base = base;

    // líquido y pared (sección rayada)
    let back = "";
    if (!only) {
      back += part("medium", rect(xl, g.yL, 2 * g.Ri, g.Hl, "bp-liq"));
      back += `<path class="bp-o" d="M${xl} ${g.yL}H${xr}" style="stroke-dasharray:none"/>`;
      const tri = [[cx + g.Ri * 0.62, g.yL], [cx + g.Ri * 0.62 - 6, g.yL - 9], [cx + g.Ri * 0.62 + 6, g.yL - 9]];
      back += `<path class="bp-o bp-f" d="${pth(tri, true)}"/>`;
      let walls = rect(xl - g.t, g.yT, g.t, g.Hh, "bp-o bp-cut") + rect(xr, g.yT, g.t, g.Hh, "bp-o bp-cut");
      back += part("shell", walls + hit(`M${xl - g.t / 2} ${g.yT + 4}V${g.yB - 4}`, 16) + hit(`M${xr + g.t / 2} ${g.yT + 4}V${g.yB - 4}`, 16));
      // deflectores (corte)
      let bf = "";
      for (const sx of [-1, 1]) bf += rect(sx < 0 ? xl : xr - g.Wb, g.yB - g.Hl - 0.02 * g.Ri, g.Wb, g.Hl + 0.02 * g.Ri - 0.025 * g.Ri * 2, "bp-o bp-cut").replace('class="bp-o bp-cut"', 'class="bp-o bp-cut bp-thin"');
      bf = bf.replace(/width="[\d.]+"/g, `width="${f1(Math.max(2.4, 0.012 * par.Di * g.s))}"`);
      // deflector: placa fina junto a la pared vista de canto, con su ancho radial como cota
      let bfs = "";
      for (const sx of [-1, 1]) {
        const w = Math.max(2.4, 0.012 * par.Di * g.s), x = sx < 0 ? xl : xr - w;
        bfs += rect(x, g.yB - g.Hl - 6, w, g.Hl + 6 - 0.025 * par.Di * g.s, "bp-o bp-f");
        bfs += `<path class="bp-thinline" d="M${f1(sx < 0 ? xl + g.Wb : xr - g.Wb)} ${f1(g.yB - g.Hl * 0.85)}H${f1(sx < 0 ? xl : xr)}" />`;
      }
      back += part("baffle", bfs + hit(`M${xl + 3} ${g.yB - 20}V${g.yL}`, 14) + hit(`M${xr - 3} ${g.yB - 20}V${g.yL}`, 14));
    }
    out.back = back;

    // serpentín: hélice vista de lado. x = Rc·cos(a), profundidad = sen(a)
    const thIn = -150, thOut = 30;
    const th0 = thOut - 360 * g.N;
    const pts = [];
    const nh = Math.max(8, Math.ceil((g.N * 360) / 5));
    for (let i = 0; i <= nh; i++) { const sN = (g.N * i) / nh, a = th0 + 360 * sN; pts.push({ x: cx + g.Rc * Math.cos(a * D2R), y: g.yB - (g.z0 + g.pv * sN), d: Math.sin(a * D2R), sN }); }
    const runs = [];
    let cur = { front: pts[0].d >= 0, pts: [pts[0]] };
    for (let i = 1; i < pts.length; i++) {
      const fr = pts[i].d >= 0;
      if (fr !== cur.front) { cur.pts.push(pts[i]); runs.push(cur); cur = { front: fr, pts: [pts[i - 1], pts[i]] }; } else cur.pts.push(pts[i]);
    }
    runs.push(cur);
    const tubeW = Math.max(2.4, g.de);
    const mk = (run) => {
      const d = pth(run.pts.map((p) => [p.x, p.y]));
      const mid = run.pts[Math.floor(run.pts.length / 2)].sN / g.N;
      const col = cScale(Tin + (Tout - Tin) * mid);
      return part("coil", `<path class="bp-glow" d="${d}" style="stroke:${col};stroke-width:${f1(tubeW * 4.2)}"/><path class="bp-tube-edge" d="${d}" style="stroke-width:${f1(tubeW + 1.8)}"/><path class="bp-tube" d="${d}" style="stroke-width:${f1(tubeW)}"/><path class="bp-core" d="${d}" style="stroke:${col};stroke-width:${f1(tubeW * 0.45)}"/><path class="bp-flowline" data-flow d="${d}" style="stroke-width:${f1(tubeW * 0.2)}"/>` + hit(d, Math.max(14, tubeW + 8)));
    };
    const backRuns = runs.filter((r) => !r.front).map(mk).join("");
    const frontRuns = runs.filter((r) => r.front).map(mk).join("");
    // tubos de entrada y salida
    const yTop = g.yT - g.lidT * 0.2;
    const xin = cx - g.Rr * Math.cos(30 * D2R), xout = cx + g.Rc * Math.cos(thOut * D2R);
    const first = pts[0], last = pts[pts.length - 1];
    const dIn = pth([[xin, yTop - 10], [xin, g.yB - g.z0], [first.x, g.yB - g.z0]]);
    const dOut = pth([[last.x, last.y], [last.x, yTop - 10]]);
    const feed = (id, d, T) => part(id, `<path class="bp-tube-edge" d="${d}" style="stroke-width:${f1(tubeW + 1.8)}"/><path class="bp-tube" d="${d}" style="stroke-width:${f1(tubeW)}"/><path class="bp-core" d="${d}" style="stroke:${cScale(T)};stroke-width:${f1(tubeW * 0.45)}"/><path class="bp-flowline" data-flow d="${d}" style="stroke-width:${f1(tubeW * 0.2)}"/>` + hit(d, Math.max(14, tubeW + 8)));
    out.coilBack = only ? "" : backRuns;
    if (only) out.coilBack = backRuns;

    // soportes (varillas)
    let rods = "";
    for (const sx of [-1, 1]) {
      const x = cx + sx * g.Rr;
      rods += `<path class="bp-o" d="M${f1(x)} ${f1(g.yT)}V${f1(g.yB - g.z0 - g.de * 1.1)}" style="stroke-width:2"/>`;
      for (let k = 1; k < g.N; k += 2) { const yy = g.yB - (g.z0 + g.pv * k); rods += `<path class="bp-o" d="M${f1(x)} ${f1(yy)}H${f1(cx + sx * (g.Rc + g.de * 0.5))}" style="stroke-width:3"/>`; }
    }
    // flujo interno: lazos
    let flow = "";
    if (!only) {
      const zc = Math.min(g.Dimp, g.Hl * 0.34);
      for (const [z, up] of [[zc, true], [zc, false]]) {
        const rl = g.Ri - g.Wb * 1.15, rc = (g.Rimp * 1.06 + rl) / 2, rx = (rl - g.Rimp * 1.06) / 2, rz = Math.min(g.Dimp * 0.42, 0.5 * (up ? g.Hl - z - 3 : z - 6));
        if (rz < 5) continue;
        for (const sx of [-1, 1]) {
          let d = "";
          for (let i = 0; i <= 36; i++) { const tt = (up ? -90 + (360 * i) / 36 : 90 - (360 * i) / 36) * D2R; d += (i ? "L" : "M") + f1(cx + sx * (rc + rx * Math.cos(tt))) + " " + f1(Y(g, z + (up ? rz : -rz) + rz * Math.sin(tt))); }
          flow += `<path class="bp-loop" data-flow d="${d}"/>`;
        }
      }
    }
    out.flow = flow;

    // eje
    const rs = Math.max(2, 0.013 * par.Di * g.s);
    out.shaft = only ? part("shaft", rect(cx - rs, g.yT - g.lidT - 14, 2 * rs, g.Hh + 14 + g.lidT - 6, "bp-o bp-f") + hit(`M${cx} ${g.yT - 10}V${g.yB - 8}`, 16)) : part("shaft", rect(cx - rs, g.yT - g.lidT - 0.2 * par.Di * g.s, 2 * rs, g.Hh + g.lidT + 0.2 * par.Di * g.s - 6, "bp-o bp-f") + hit(`M${cx} ${g.yT - 10}V${g.yB - 8}`, 16));

    // frente: varillas y serpentín delantero
    out.front = part("support", rods + hit(`M${cx - g.Rr} ${g.yT + 4}V${g.yB - g.z0}`, 12) + hit(`M${cx + g.Rr} ${g.yT + 4}V${g.yB - g.z0}`, 12)) + frontRuns + feed("cin", dIn, Tin) + feed("cout", dOut, Tout);

    // flechas de calor
    let heat = "";
    if (!only) {
      for (const [frac, side] of [[0.3, -1], [0.5, 1], [0.7, -1], [0.85, 1]]) {
        const sN = g.N * frac, y = g.yB - (g.z0 + g.pv * sN);
        const x1 = cx + side * (g.Rc + g.de / 2 + (g.Ri - g.Wb - g.Rc) * 0.45), x2 = cx + side * (g.Rc + g.de * 0.9);
        heat += `<path class="bp-heat" data-flow d="M${f1(x1)} ${f1(y - g.pv * 0.3)}L${f1(x2)} ${f1(y - 1)}"/>`;
        const xi = cx - side * (g.Rimp + (g.Rc - g.de / 2 - g.Rimp) * 0.4), xi2 = cx - side * (g.Rc - g.de * 0.9);
        heat += `<path class="bp-heat" data-flow d="M${f1(xi)} ${f1(y - g.pv * 0.3)}L${f1(xi2)} ${f1(y - 1)}"/>`;
      }
    }
    out.heat = `<g class="bp-heatg">${heat}</g>`;
    out.veil = "";

    // tapa, motor, sondas, boquillas
    let top = "";
    const Rl = g.Ro + 0.035 * par.Di * g.s;
    if (!only) {
      top += part("lid", rect(cx - Rl, g.yT - g.lidT, 2 * Rl, g.lidT, "bp-o bp-cut") + rect(cx - Rl, g.yT, 2 * Rl, 0.03 * par.Di * g.s, "bp-o bp-f") + hit(`M${cx - Rl} ${g.yT - g.lidT / 2}H${cx + Rl}`, 14));
      const Di = par.Di * g.s, yL0 = g.yT - g.lidT;
      let motor = rect(cx - 0.075 * Di, yL0 - 0.06 * Di, 0.15 * Di, 0.06 * Di, "bp-o bp-f") + rect(cx - 0.1 * Di, yL0 - 0.13 * Di, 0.2 * Di, 0.07 * Di, "bp-o bp-f") + rect(cx - 0.115 * Di, yL0 - 0.33 * Di, 0.23 * Di, 0.2 * Di, "bp-o bp-f") + rect(cx - 0.1 * Di, yL0 - 0.365 * Di, 0.2 * Di, 0.035 * Di, "bp-o bp-f");
      for (let i = -3; i <= 3; i++) motor += `<path class="bp-thinline" d="M${f1(cx + i * 0.032 * Di)} ${f1(yL0 - 0.31 * Di)}V${f1(yL0 - 0.15 * Di)}"/>`;
      top += part("motor", motor + hit(`M${cx} ${yL0 - 4}V${yL0 - 0.36 * Di}`, 40));
      // sondas
      const pr = (id, x, w, yEnd) => part(id, rect(x - w / 2, yEnd, w, g.yT - yEnd, "bp-o bp-f") + rect(x - w * 1.1, g.yT - g.lidT - 0.045 * Di, w * 2.2, 0.045 * Di + g.lidT, "bp-o bp-f") + hit(`M${x} ${yEnd}V${g.yT - g.lidT}`, 14));
      const xp1 = cx + (g.Ri - g.Wb - g.Rc * 0.0 - (g.Ri - g.Wb - (g.Rr + g.de / 2)) / 2) * -1, xp2 = cx + (g.Ri - g.Wb - (g.Ri - g.Wb - (g.Rr + g.de / 2)) / 2) * 0.998;
      top += pr("tprobe", xp1 + 0, Math.max(2, 0.009 * Di * 2), g.yB - g.Hl * 0.42) + pr("phprobe", xp2, Math.max(3.4, 0.015 * Di * 2), g.yB - g.Hl * 0.5);
      // boquillas del medio
      const nz = (id, y, dir) => { const x0 = dir < 0 ? xl - g.t : xr + g.t, len = 0.16 * Di, x1 = x0 + dir * len, r = 0.04 * Di; return part(id, rect(Math.min(x0, x1), y - r, len, 2 * r, "bp-o bp-f") + rect(x1 - (dir > 0 ? 0 : 0.018 * Di), y - r * 1.7, 0.018 * Di, r * 3.4, "bp-o bp-f") + `<path class="bp-medline" data-flow d="M${f1(x1 + dir * (dir > 0 ? 0.2 : 0.02) * Di * (dir > 0 ? 1 : -1) * (dir > 0 ? 1 : -1))} ${y}H${f1(x1 + dir * 0.02 * Di)}"/>` + hit(`M${f1(x0)} ${y}H${f1(x1)}`, 18)); };
      top += nz("inlet", g.yB - g.Hh * 0.8, -1) + nz("outlet", g.yB - (g.Hh * 0.12 + 6), 1);
    }
    // conexiones del refrigerante sobre la tapa
    const yF = g.yT - g.lidT;
    const fit = (id, x, T, dir) => { const h = g.de * 2.6, r = g.de * 0.62; return part(id, rect(x - g.de * 1.5, yF - 2, g.de * 3, 2, "bp-o bp-f") + rect(x - r, yF - 2 - h, 2 * r, h, "bp-o bp-f") + `<path class="bp-tube-edge" d="M${f1(x)} ${f1(yF - 2 - h)}H${f1(x + dir * 0.1 * par.Di * g.s)}" style="stroke-width:${f1(tubeW + 1.8)}"/><path class="bp-tube" d="M${f1(x)} ${f1(yF - 2 - h)}H${f1(x + dir * 0.1 * par.Di * g.s)}" style="stroke-width:${f1(tubeW)}"/><path class="bp-core" d="M${f1(x)} ${f1(yF - 2 - h)}H${f1(x + dir * 0.1 * par.Di * g.s)}" style="stroke:${cScale(T)};stroke-width:${f1(tubeW * 0.45)}"/>` + hit(`M${f1(x)} ${f1(yF - 2)}V${f1(yF - 2 - h)}`, 18)); };
    top += fit("cin", xin, Tin, -1) + fit("cout", xout, Tout, 1);
    out.shell = "";

    // cotas
    let dim = "";
    if (dims) {
      const arrow = (x1, y1, x2, y2) => { const a = Math.atan2(y2 - y1, x2 - x1), h = 4.5; const hd = (x, y, s) => `M${f1(x + s * h * Math.cos(a - 0.4))} ${f1(y + s * h * Math.sin(a - 0.4))}L${f1(x)} ${f1(y)}L${f1(x + s * h * Math.cos(a + 0.4))} ${f1(y + s * h * Math.sin(a + 0.4))}`; return `<path class="bp-dimline" d="M${f1(x1)} ${f1(y1)}L${f1(x2)} ${f1(y2)}${hd(x1, y1, 1)}${hd(x2, y2, -1)}"/>`; };
      const txt = (x, y, t, anc = "middle", rot = 0) => `<text class="bp-dimtext" x="${f1(x)}" y="${f1(y)}" text-anchor="${anc}"${rot ? ` transform="rotate(${rot} ${f1(x)} ${f1(y)})"` : ""}>${t}</text>`;
      const nf1 = (v) => v.toLocaleString("es-MX", { maximumFractionDigits: 1 });
      // Di
      const yD = g.yB + g.tb + g.leg * 0.45;
      dim += `<path class="bp-ext" d="M${xl} ${g.yB + 2}V${yD + 6}M${xr} ${g.yB + 2}V${yD + 6}"/>` + arrow(xl, yD, xr, yD) + txt(cx, yD - 4, `Di ${nf1(par.Di)} cm`);
      // Ht
      const xH = xl - g.Ro * 0.62;
      dim += `<path class="bp-ext" d="M${xl - g.t - 3} ${g.yT}H${xH - 6}M${xl - g.t - 3} ${g.yB}H${xH - 6}"/>` + arrow(xH, g.yT, xH, g.yB) + txt(xH - 4, (g.yT + g.yB) / 2, `Ht ${nf1(par.Ht)} cm`, "middle", -90);
      // Hl
      const xHl = xr + g.Ro * 0.62;
      dim += `<path class="bp-ext" d="M${xr + g.t + 3} ${g.yL}H${xHl + 6}M${xr + g.t + 3} ${g.yB}H${xHl + 6}"/>` + arrow(xHl, g.yL, xHl, g.yB) + txt(xHl + 12, (g.yL + g.yB) / 2, `Hl ${nf1(g.Hl / g.s)} cm`, "middle", 90);
      // Dc
      const yDc = g.yB - g.z0 + g.de * 2.6;
      dim += `<path class="bp-ext" d="M${cx - g.Rc} ${g.yB - g.z0}V${yDc + 4}M${cx + g.Rc} ${g.yB - g.z0}V${yDc + 4}"/>` + arrow(cx - g.Rc, yDc, cx + g.Rc, yDc) + txt(cx, yDc - 4, `Dc ${nf1(par.Dc)} cm`);
      // paso p entre dos espiras a la derecha
      const yp1 = g.yB - (g.z0 + g.pv * (g.N * 0.5)), yp2 = yp1 - g.pv, xpp = cx + g.Rc + g.de * 3.2;
      dim += `<path class="bp-ext" d="M${cx + g.Rc + g.de} ${yp1}H${xpp + 5}M${cx + g.Rc + g.de} ${yp2}H${xpp + 5}"/>` + arrow(xpp, yp2, xpp, yp1) + txt(xpp + 8, yp1 - g.pv / 2 + 3, `p ${nf1(par.p)} cm`, "start");
      // altura de la hélice N·p
      const yc0 = g.yB - g.z0, yc1 = g.yB - (g.z0 + g.pv * g.N), xC = cx - g.Rc - g.de * 3.4;
      dim += `<path class="bp-ext" d="M${cx - g.Rc - g.de} ${yc0}H${xC - 4}M${cx - g.Rc - g.de} ${yc1}H${xC - 4}"/>` + arrow(xC, yc1, xC, yc0) + txt(xC - 4, (yc0 + yc1) / 2, `N·p ${nf1(par.N * par.p)} cm`, "middle", -90);
    }
    out.top = top + (dims ? `<g class="bp-dims">${dim}</g>` : "");
    out.imp = { g, par };

    // anclas de etiquetas
    anchors.coil = { p: [cx + g.Rc, g.yB - (g.z0 + g.pv * g.N * 0.4)], side: "R", text: "Serpentín helicoidal" };
    anchors.cin = { p: [xin, yF - 2 - g.de * 2.6], side: "L", text: "Entrada refrigerante" };
    anchors.cout = { p: [xout, yF - 2 - g.de * 2.6], side: "R", text: "Salida refrigerante" };
    anchors.support = { p: [cx + g.Rr, g.yT + g.Hh * 0.12], side: "R", text: "Soporte" };
    anchors.shaft = { p: [cx, g.yT + g.Hh * 0.08], side: "L", text: "Eje de agitación" };
    anchors.impeller = { p: [cx + g.Rimp * 0.8, g.yB - Math.min(g.Dimp, g.Hl * 0.34)], side: "L", text: "Impulsor" };
    anchors.baffle = { p: [xl + 3, g.yB - g.Hl * 0.6], side: "L", text: "Deflector" };
    anchors.medium = { p: [cx - g.Ri * 0.6, g.yL + 8], side: "L", text: "Medio de cultivo" };
    anchors.shell = { p: [xr + g.t / 2, g.yT + g.Hh * 0.5], side: "R", text: "Pared del tanque" };
    anchors.lid = { p: [cx - Rl * 0.8, g.yT - g.lidT / 2], side: "L", text: "Tapa" };
    anchors.motor = { p: [cx + 0.115 * par.Di * g.s, g.yT - g.lidT - 0.2 * par.Di * g.s], side: "R", text: "Motor y reductor" };
    anchors.inlet = { p: [xl - g.t - 0.14 * par.Di * g.s, g.yB - g.Hh * 0.8], side: "L", text: "Entrada del medio" };
    anchors.outlet = { p: [xr + g.t + 0.14 * par.Di * g.s, g.yB - (g.Hh * 0.12 + 6)], side: "R", text: "Salida del medio" };
    return { g, out, anchors, par };
  }

  // Impulsores en vista lateral: las palas se ven con ancho r·cos(ángulo)
  function dynamic(b, tsec) {
    const g = b.out.imp.g, cx = g.cx;
    const ang = tsec * 360 * 0.42;
    const zc1 = Math.min(g.Dimp, g.Hl * 0.34), zc2 = zc1 + 1.5 * g.Dimp;
    const imps = [zc1].concat(zc2 + g.Dimp * 0.24 + 4 < g.Hl ? [zc2] : []);
    let s = "";
    for (const zc of imps) {
      const y = g.yB - zc, Rd = g.Rimp * 0.74, h = Math.max(4, g.Rimp * 0.8), r0 = g.Rimp * 0.48;
      s += rect(cx - Rd, y - 1.4, 2 * Rd, 2.8, "bp-o bp-f");
      for (let k = 0; k < 6; k++) {
        const a = (ang + k * 60) * D2R, c = Math.cos(a), sn = Math.sin(a);
        const xa = cx + r0 * c, xb = cx + g.Rimp * c;
        const th = Math.max(1.2, 0.012 * g.Rimp * 2 * Math.abs(sn) + 1.2);
        s += `<path class="bp-blade${sn < 0 ? " bp-back" : ""}" d="M${f1(xa)} ${f1(y - h / 2)}H${f1(xb)}V${f1(y + h / 2)}H${f1(xa)}Z" style="opacity:${sn < 0 ? 0.5 : 1}"/>`;
      }
      s += rect(cx - 2.2, y - h / 2, 4.4, h, "bp-o bp-f");
    }
    return { back: "", mid: part("impeller", s + hit(`M${cx} ${g.yB - zc1 - 14}V${g.yB - zc1 + 14}`, 30)), front: "", ang };
  }

  function labels(b, fontPx, ids) {
    const fs = (fontPx * W) / (document.querySelector("#fig svg")?.getBoundingClientRect().width || W);
    const list = ids.map((id) => ({ id, ...b.anchors[id] })).filter((q) => q.p);
    for (const q of list) { q.s = q.p; q.side = q.s[0] < W / 2 ? "L" : "R"; }
    const place = (side) => {
      const mine = list.filter((q) => q.side === side).sort((a, c) => a.s[1] - c.s[1]);
      const gap = fs * 1.9;
      let y = -1e9;
      for (const q of mine) { q.y = Math.max(q.s[1], y + gap); y = q.y; }
      const over = y - (H - 14);
      if (over > 0) for (const q of mine) q.y -= over;
      let low = 12 + fs;
      for (const q of mine) { q.y = Math.max(q.y, low); low = q.y + gap; }
    };
    place("L"); place("R");
    let s = "";
    for (const q of list) {
      const left = q.side === "L", x = left ? 8 : W - 8;
      const dxAbs = Math.abs(q.y - q.s[1]);
      const mx = left ? Math.max(x + 4, q.s[0] - dxAbs) : Math.min(x - 4, q.s[0] + dxAbs);
      const end = left ? x + q.text.length * fs * 0.56 + 4 : x - q.text.length * fs * 0.56 - 4;
      s += `<g class="bx-label" data-label="${q.id}"><path class="bx-leader" d="M${f1(q.s[0])} ${f1(q.s[1])}L${f1(mx)} ${f1(q.y)}L${f1(end)} ${f1(q.y)}"/><circle class="bx-leader-dot" cx="${f1(q.s[0])}" cy="${f1(q.s[1])}" r="1.9"/><text x="${x}" y="${f1(q.y + fs * 0.34)}" font-size="${fs.toFixed(1)}" text-anchor="${left ? "start" : "end"}">${q.text}</text></g>`;
    }
    return s;
  }
  return { W, get H() { return H; }, build, dynamic, labels, tempColor };
})();

// ===== Interfaz: estado, paneles, cálculos y animación =====
const APP = (() => {
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const nf = (n, d = 2) => (Number.isFinite(n) ? n.toLocaleString("es-MX", { minimumFractionDigits: d, maximumFractionDigits: d }) : "—");
  const clampn = (v, a, b) => Math.min(b, Math.max(a, v));

  const VAR = typeof VARIANT !== "undefined" ? VARIANT : {};
  const DEF = { ...VALORES_ILUSTRATIVOS, ...(VAR.defaults || {}) };
  const state = {
    par: sincronizar({ ...DEF, L: 0 }, "N"),
    anchor: "N",
    view: "cut",
    labels: true,
    tab: "comp",
    sel: null,
    tshow: DEF.T0,
    paused: false,
    built: null,
    goodPar: null,
    invalid: false,
  };
  state.goodPar = { ...state.par };
  const still = matchMedia("(prefers-reduced-motion: reduce)");
  const G = { dyn: {}, t0: performance.now(), last: 0, visible: true, raf: 0 };

  // ---------- Dibujo ----------
  function layers() {
    return { base: $("#g-base"), back: $("#g-back"), coilBack: $("#g-coilBack"), flow: $("#g-flow"), dynBack: $("#g-dynBack"), shaft: $("#g-shaft"), dyn: $("#g-dyn"), dynFront: $("#g-dynFront"), front: $("#g-front"), heat: $("#g-heat"), veil: $("#g-veil"), shell: $("#g-shell"), top: $("#g-top"), labels: $("#g-labels") };
  }
  function rebuild() {
    const par = state.goodPar;
    const compact = ($("#svg").getBoundingClientRect().width || 600) < 520 && !state.labels;
    state.compact = compact;
    const b = SCENE.build(par, { view: state.view, compact });
    $("#svg").setAttribute("viewBox", `0 0 ${SCENE.W} ${SCENE.H}`);
    state.built = b;
    const l = layers();
    l.base.innerHTML = b.out.base;
    l.back.innerHTML = b.out.back;
    l.coilBack.innerHTML = b.out.coilBack;
    l.flow.innerHTML = b.out.flow;
    l.shaft.innerHTML = b.out.shaft;
    l.front.innerHTML = b.out.front;
    l.heat.innerHTML = b.out.heat;
    l.veil.innerHTML = b.out.veil;
    l.shell.innerHTML = b.out.shell;
    l.top.innerHTML = b.out.top;
    paintTemps();
    frame(performance.now(), true);
    drawLabels();
    markSelection();
    readout();
  }
  function paintTemps() {
    const p = state.goodPar;
    const svg = $("#svg");
    const u = (T) => (T - p.Tin) / Math.max(1, p.T0 - p.Tin);
    const c = SCENE.tempColor(u(state.tshow));
    svg.style.setProperty("--bx-medium", c);
    const dur = clampn(2.4 * Math.sqrt(11 / Math.max(0.5, p.caudal)), 0.7, 3.6);
    svg.style.setProperty("--bx-flow-dur", dur.toFixed(2) + "s");
  }
  function drawLabels() {
    if (!state.built) return;
    const l = layers();
    if (!state.labels) { l.labels.innerHTML = ""; return; }
    const w = $("#svg").getBoundingClientRect().width || SCENE.W;
    const all = ["lid", "shaft", "impeller", "tprobe", "phprobe", "outlet", "cin", "medium", "coil", "cout", "support", "baffle", "motor", "inlet", "shell", "drain"];
    const small = ["cin", "impeller", "tprobe", "medium", "coil", "cout", "motor", "drain"];
    const ids = VAR.labels || (w < 520 ? small : all);
    l.labels.innerHTML = SCENE.labels(state.built, w < 520 ? 11 : 11.5, ids);
    markSelection();
  }
  function frame(now, force) {
    if (!state.built) return;
    const t = (now - G.t0) / 1000;
    const d = SCENE.dynamic(state.built, t);
    const l = layers();
    l.dynBack.innerHTML = d.back;
    l.dyn.innerHTML = d.mid;
    l.dynFront.innerHTML = d.front;
    if (state.sel) markSelection();
  }
  function loop(now) {
    G.raf = 0;
    if (state.paused || still.matches || !G.visible || document.hidden) return;
    if (now - G.last >= 1000 / 30) { G.last = now; frame(now); }
    G.raf = requestAnimationFrame(loop);
  }
  function kick() { if (!G.raf && !state.paused && !still.matches && G.visible && !document.hidden) G.raf = requestAnimationFrame(loop); }

  // ---------- Selección ----------
  function select(id) {
    state.sel = id && COMPONENTES[id] ? id : null;
    markSelection();
    renderInfo();
    if (state.sel && state.tab !== "comp") setTab("comp");
  }
  function markSelection() {
    $$("[data-part].sel,.bx-label.sel").forEach((e) => e.classList.remove("sel"));
    if (state.sel) $$(`[data-part="${state.sel}"],.bx-label[data-label="${state.sel}"]`).forEach((e) => e.classList.add("sel"));
    $$(".chip").forEach((c) => c.setAttribute("aria-pressed", c.dataset.id === state.sel ? "true" : "false"));
    const r = $("[data-readout]");
    if (r && state.sel) r.textContent = COMPONENTES[state.sel].nombre;
    else readout();
  }
  function readout() {
    const r = $("[data-readout]");
    if (!r) return;
    const c = calcular(state.goodPar, state.tshow);
    r.textContent = c.ok ? `L ${nf(c.L, 2)} m · A ${nf(c.A, 3)} m² · Q̇ ${nf(c.act.q, 0)} W` : "revisa los parámetros";
  }
  function paramRow(key) {
    const c = CAMPOS[key];
    const v = state.par[key];
    return `<li><span>${c.etiqueta}</span><b>${nf(v, key === "N" ? 2 : 1)} ${c.unidad}</b></li>`;
  }
  function renderInfo() {
    const box = $("#info");
    const id = state.sel;
    if (!id) {
      box.innerHTML = `<p class="hint">Toca un componente en la ilustración o elígelo de la lista para ver su función, material y parámetros.</p>`;
      return;
    }
    const c = COMPONENTES[id];
    const lista = c.params.length ? `<ul class="plist">${c.params.map(paramRow).join("")}</ul>` : `<p class="muted">Este componente no depende de los parámetros del panel.</p>`;
    box.innerHTML = `<h3>${c.nombre}</h3>
      <dl class="facts">
        <dt>Función</dt><dd>${c.funcion}</dd>
        <dt>Descripción técnica</dt><dd>${c.descripcion}</dd>
        <dt>Material sugerido</dt><dd>${c.material}</dd>
        <dt>Importancia en el sistema</dt><dd>${c.importancia}</dd>
      </dl>
      <h4>Parámetros relacionados</h4>${lista}`;
  }

  // ---------- Parámetros ----------
  const GRUPOS = [
    { id: "geo", titulo: "Geometría", claves: ["Vw", "Di", "Ht", "De", "Dc", "N", "L", "p"] },
    { id: "term", titulo: "Condiciones térmicas", claves: ["T0", "Tobj", "Tin", "Tout", "caudal", "U"] },
  ];
  function buildForm() {
    $("#form").innerHTML = GRUPOS.map((g) => `<fieldset><legend>${g.titulo}</legend>${g.claves.map((k) => {
      const c = CAMPOS[k];
      const note = k === "L" ? `<span class="note">Se calcula con L = N·√[(π·Dc)² + p²]. También puedes escribirla y se ajustan las vueltas.</span>` : k === "Vw" ? `<span class="note">Volumen de trabajo; el nivel de líquido sale de Vw y Di.</span>` : "";
      return `<div class="field" data-k="${k}"><label for="n-${k}">${c.etiqueta}</label><div class="inp"><input id="n-${k}" type="text" inputmode="decimal" autocomplete="off" value="${fmtIn(state.par[k], k)}" data-k="${k}"><span class="u">${c.unidad}</span></div><input class="rng" type="range" min="${c.min}" max="${c.max}" step="${c.paso}" value="${clampn(state.par[k], c.min, c.max)}" data-r="${k}" aria-label="${c.etiqueta}">${note}<div class="err" id="e-${k}" role="alert"></div></div>`;
    }).join("")}</fieldset>`).join("");
    $$("#form input[type=text]").forEach((i) => i.addEventListener("input", () => onType(i.dataset.k, i.value)));
    $$("#form input[type=range]").forEach((i) => i.addEventListener("input", () => { const k = i.dataset.r; $(`#n-${k}`).value = fmtIn(Number(i.value), k); onType(k, i.value, true); }));
  }
  const fmtIn = (v, k) => (Number.isFinite(v) ? String(Math.round(v * (k === "N" ? 100 : 1000)) / (k === "N" ? 100 : 1000)) : "");
  function onType(k, raw, fromRange) {
    const v = Number(String(raw).trim().replace(",", "."));
    const par = { ...state.par, [k]: raw === "" ? NaN : v };
    if (k === "L") state.anchor = "L";
    else if (k === "N") state.anchor = "N";
    if (Number.isFinite(par.Dc) && Number.isFinite(par.p) && par.Dc > 0 && par.p > 0) {
      if (k === "L" || (state.anchor === "L" && k !== "N")) { if (Number.isFinite(par.L)) par.N = sincronizar(par, "L").N; }
      else if (Number.isFinite(par.N)) par.L = sincronizar(par, "N").L;
    }
    state.par = par;
    for (const o of ["N", "L"]) if (o !== k) { const inp = $(`#n-${o}`); if (inp && Number.isFinite(par[o])) inp.value = fmtIn(par[o], o); const r = $(`[data-r="${o}"]`); if (r) r.value = clampn(par[o], CAMPOS[o].min, CAMPOS[o].max); }
    if (!fromRange) { const r = $(`[data-r="${k}"]`); if (r && Number.isFinite(v)) r.value = clampn(v, CAMPOS[k].min, CAMPOS[k].max); }
    schedule();
  }
  let pend = 0;
  function schedule() { if (!pend) pend = requestAnimationFrame(() => { pend = 0; applyParams(); }); }
  function applyParams() {
    const v = validar(state.par);
    $$(".err").forEach((e) => (e.textContent = ""));
    $$(".field").forEach((f) => f.classList.remove("bad", "warn"));
    const mark = (list, cls) => list.forEach((e) => { const el = $(`#e-${e.campo}`); if (el && !el.textContent) { el.textContent = e.msg; el.dataset.kind = cls; } $(`.field[data-k="${e.campo}"]`)?.classList.add(cls); });
    mark(v.avisos, "warn");
    mark(v.errores, "bad");
    const banner = $("#banner");
    if (v.errores.length) {
      banner.hidden = false;
      banner.textContent = `Hay ${v.errores.length} ${v.errores.length === 1 ? "valor" : "valores"} por corregir. La ilustración conserva la última configuración válida y los cálculos esperan.`;
      state.invalid = true;
    } else {
      banner.hidden = true;
      state.invalid = false;
      const prev = state.goodPar;
      state.goodPar = { ...state.par };
      const Tmin = state.par.Tobj, Tmax = state.par.T0;
      state.tshow = clampn(state.tshow, Tmin, Tmax);
      if (prev.T0 !== state.par.T0) state.tshow = state.par.T0;
      syncT();
      rebuild();
    }
    renderCalc();
    emit();
  }
  function syncT() {
    const r = $("#tshow");
    r.min = state.goodPar.Tobj; r.max = state.goodPar.T0; r.step = 0.5; r.value = state.tshow;
    $("#tshowv").textContent = nf(state.tshow, 1) + " °C";
  }

  // ---------- Cálculos ----------
  const row = (k, f, v) => `<div class="row"><span class="k">${k}</span><code class="f">${f}</code><span class="v">${v}</span></div>`;
  function renderCalc() {
    const box = $("#calc");
    if (state.invalid) {
      const v = validar(state.par);
      box.innerHTML = `<p class="hint">Los cálculos se muestran cuando todos los valores son válidos.</p><ul class="errs">${v.errores.map((e) => `<li>${e.msg}</li>`).join("")}</ul>`;
      return;
    }
    const c = calcular(state.goodPar, state.tshow);
    const p = state.goodPar;
    const g = c.geo;
    const estado = { coherente: "Los dos balances coinciden dentro de la tolerancia.", revisar: "Los balances difieren; conviene ajustar caudal, T de salida, área o U.", inconsistente: "Los balances no son coherentes entre sí: revisa los datos de entrada." }[c.estado];
    box.innerHTML = `
      <section><h4>Geometría del serpentín</h4>
        ${row("Longitud de la hélice", "L = N·√[(π·Dc)² + p²]", `${nf(c.L, 2)} m`)}
        ${row("Superficie de transferencia", "A = π·De·L", `${nf(c.A, 3)} m²`)}
        ${row("Altura ocupada por la hélice", "N·p", `${nf(g.altura * 100, 1)} cm`)}
        ${row("Nivel de líquido", "Hl = 4·Vw/(π·Di²)", `${nf(g.Hl * 100, 1)} cm (${nf(c.llenado * 100, 0)} % de la altura)`)}
      </section>
      <section><h4>Calor sensible del cultivo</h4>
        ${row("Masa del medio", "m = ρ·Vw", `${nf(c.masa, 1)} kg`)}
        ${row("Energía a retirar", "Q = m·Cp·(T0 − Tobj)", `${nf(c.energia / 1000, 0)} kJ`)}
      </section>
      <section><h4>Transferencia por el serpentín (T del cultivo: ${nf(c.T, 1)} °C)</h4>
        ${row("Diferencias en los extremos", "ΔT1 = Tc − Tent, ΔT2 = Tc − Tsal", `${nf(c.act.dT1, 1)} K y ${nf(c.act.dT2, 1)} K`)}
        ${row("Media logarítmica", "ΔTlm = (ΔT1 − ΔT2)/ln(ΔT1/ΔT2)", `${nf(c.act.dtlm, 2)} K`)}
        ${row("Potencia transferible", "Q̇ = U·A·ΔTlm", `${nf(c.act.q, 0)} W`)}
        ${row("Con T inicial / con T objetivo", "mismo cálculo con Tc = T0 y Tc = Tobj", `${nf(c.ini.q, 0)} W / ${nf(c.obj.q, 0)} W`)}
      </section>
      <section><h4>Balance del refrigerante</h4>
        ${row("Gasto másico", "ṁ = caudal·ρ", `${nf(c.gasto, 3)} kg/s`)}
        ${row("Potencia absorbida", "q = ṁ·Cp·(Tsal − Tent)", `${nf(c.qRef, 0)} W`)}
        ${row("Diferencia entre balances", "|q − Q̇|/Q̇", `${nf(c.dif * 100, 1)} %`)}
        ${row("Área que haría falta", "Areq = q/(U·ΔTlm)", `${nf(c.areaNecesaria, 3)} m² (hay ${nf(c.A, 3)} m²)`)}
        <p class="status" data-s="${c.estado}">${estado}</p>
      </section>
      <details open><summary>Supuestos</summary><ul class="sup">
        <li>El medio y el refrigerante se tratan como agua: ρ = ${SUPUESTOS.rhoMedio} kg/m³ y Cp = ${SUPUESTOS.cpMedio} J/(kg·K) para el medio; Cp = ${SUPUESTOS.cpRef} J/(kg·K) para el refrigerante.</li>
        <li>El cultivo está bien mezclado y a una sola temperatura Tc. El refrigerante entra a Tent y sale a Tsal.</li>
        <li>U es un dato que escribes tú. No se calcula con correlaciones y se supone constante a lo largo de la tubería.</li>
        <li>A usa solo la superficie exterior de la hélice. No cuenta los tramos rectos de entrada y salida ni los soportes.</li>
        <li>No hay generación de calor metabólico, pérdidas al ambiente ni calor de agitación. El resultado es un cálculo de dimensionamiento aproximado y no una simulación del proceso.</li>
        <li>Las proporciones del impulsor (Di/3) y de los deflectores (Di/10) son referencias de dibujo, no datos de diseño.</li>
        <li>Los valores iniciales de la pantalla son ilustrativos y no corresponden a ningún experimento.</li>
      </ul></details>`;
  }

  // ---------- Pestañas, chips y controles ----------
  function setTab(t) {
    state.tab = t;
    $$(".tab").forEach((b) => b.setAttribute("aria-selected", b.dataset.tab === t ? "true" : "false"));
    $$(".panel").forEach((p) => (p.hidden = p.dataset.tab !== t));
  }
  function emit() {
    const r = calcular(state.goodPar, state.tshow);
    const detail = { parametros: { ...state.goodPar }, valido: !state.invalid, resultados: r.ok ? { L_m: r.L, A_m2: r.A, energia_J: r.energia, Q_W: r.act.q, q_refrigerante_W: r.qRef, estado: r.estado } : null };
    try { window.dispatchEvent(new CustomEvent("biorreactor:cambio", { detail })); } catch (e) {}
    try { if (window.AndroidBridge && typeof window.AndroidBridge.onCambio === "function") window.AndroidBridge.onCambio(JSON.stringify(detail)); } catch (e) {}
  }
  function init() {
    $("#chips").innerHTML = ORDEN_COMPONENTES.map((id) => `<button class="chip" type="button" data-id="${id}" aria-pressed="false">${COMPONENTES[id].nombre}</button>`).join("");
    $("#chips").addEventListener("click", (e) => { const b = e.target.closest(".chip"); if (b) select(state.sel === b.dataset.id ? null : b.dataset.id); });
    $$(".tab").forEach((b) => b.addEventListener("click", () => setTab(b.dataset.tab)));
    $$("[data-view]").forEach((b) => b.addEventListener("click", () => { state.view = b.dataset.view; $$("[data-view]").forEach((x) => x.setAttribute("aria-pressed", x === b ? "true" : "false")); rebuild(); }));
    $("#btnLabels").addEventListener("click", (e) => { state.labels = !state.labels; e.currentTarget.setAttribute("aria-pressed", String(state.labels)); rebuild(); });
    $("#btnPause").addEventListener("click", (e) => { state.paused = !state.paused; e.currentTarget.setAttribute("aria-pressed", String(state.paused)); e.currentTarget.textContent = state.paused ? "Reanudar" : "Pausar"; $("#svg").classList.toggle("paused", state.paused); kick(); });
    $("#btnTheme").addEventListener("click", () => { const light = document.body.dataset.theme !== "light"; document.body.dataset.theme = light ? "light" : "dark"; });
    $("#btnReset").addEventListener("click", () => { state.par = sincronizar({ ...DEF, L: 0 }, "N"); state.anchor = "N"; state.tshow = state.par.T0; buildForm(); applyParams(); });
    $("#tshow").addEventListener("input", (e) => { state.tshow = Number(e.target.value); $("#tshowv").textContent = nf(state.tshow, 1) + " °C"; paintTemps(); renderCalc(); readout(); });
    $("#svg").addEventListener("click", (e) => { const part = e.target.closest("[data-part]"); const lab = e.target.closest("[data-label]"); const id = part?.dataset.part || lab?.dataset.label || null; select(id === state.sel ? null : id); });
    new ResizeObserver(() => { clearTimeout(init.rz); init.rz = setTimeout(() => { const small = ($("#svg").getBoundingClientRect().width || 600) < 520 && !state.labels; if (small !== state.compact) rebuild(); else drawLabels(); }, 80); }).observe($("#svg"));
    new IntersectionObserver((en) => { G.visible = en[en.length - 1].isIntersecting; if (G.visible) kick(); }, { rootMargin: "100px" }).observe($("#svg"));
    document.addEventListener("visibilitychange", kick);
    still.addEventListener?.("change", () => { if (still.matches) frame(performance.now()); kick(); });
    document.body.dataset.theme = VAR.theme || (matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark");
    if (VAR.vistas) $$("[data-view]").forEach((b) => (b.textContent = VAR.vistas[b.dataset.view] || b.textContent));
    if (($("#svg").getBoundingClientRect().width || 600) < 520) { state.labels = false; $("#btnLabels").setAttribute("aria-pressed", "false"); }
    buildForm();
    syncT();
    setTab("comp");
    applyParams();
    renderInfo();
    kick();
    window.BiorreactorApp = {
      getParametros: () => ({ ...state.goodPar }),
      setParametros: (obj) => { state.par = sincronizar({ ...state.par, ...obj }, obj && "L" in obj && !("N" in obj) ? "L" : "N"); buildForm(); applyParams(); },
      seleccionar: select,
      resultados: () => calcular(state.goodPar, state.tshow),
    };
  }
  return { init, state };
})();
APP.init();
</script></body></html>
)PLANO";

const std::string PAGINA_HTML = R"HTML(
<!doctype html>
<html lang="es">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>NuCoil | Ingeniería de Serpentines en Biorreactores</title>
  <style>
    :root {
      --bg: #f4f1e8;
      --panel: #ebe6d6;
      --panel-alt: #19243e;
      --border: #263554;
      --text: #16201f;
      --muted: #5a6b69;
      --cyan: #1f6b6a;
      --amber: #b08440;
      --emerald: #1f6b6a;
      --rose: #f43f5e;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: var(--bg);
      color: var(--text);
      line-height: 1.45;
    }
    header {
      background: #fbfaf4;
      border-bottom: 1px solid var(--border);
      padding: 14px 24px;
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 12px;
      position: sticky;
      top: 0;
      z-index: 10;
    }
    .brand { display: flex; align-items: center; gap: 12px; }
    .badge-cpp {
      background: linear-gradient(135deg, #1f6b6a, #2563eb);
      color: #fff;
      font-weight: 800;
      padding: 6px 10px;
      border-radius: 8px;
      font-size: 0.85rem;
      letter-spacing: 0.5px;
    }
    .brand h1 { margin: 0; font-size: 1.1rem; }
    .brand small { color: var(--muted); display: block; font-size: 0.78rem; }
    .toolbar { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
    button {
      background: var(--panel-alt);
      color: var(--text);
      border: 1px solid var(--border);
      padding: 7px 12px;
      border-radius: 8px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      transition: 0.15s;
    }
    button:hover { border-color: var(--cyan); color: #fff; }
    button.primary { background: #1f6b6a; border-color: #1f6b6a; color: #fff; }
    button.primary:hover { background: #0369a1; }
    .container {
      max-width: 1420px;
      margin: 0 auto;
      padding: 20px;
      display: grid;
      grid-template-columns: 440px 1fr;
      gap: 20px;
    }
    @media (max-width: 1020px) { .container { grid-template-columns: 1fr; } }
    .card {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 18px;
    }
    .card h2 {
      margin: 0 0 12px 0;
      font-size: 0.95rem;
      text-transform: uppercase;
      letter-spacing: 0.6px;
      color: var(--cyan);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    .tabs { display: flex; gap: 6px; margin-bottom: 14px; }
    .tab-btn { flex: 1; text-align: center; padding: 8px 6px; font-size: 0.78rem; }
    .tab-btn.active { background: #1f6b6a; border-color: #1f6b6a; color: #fff; }
    .tab-pane { display: none; }
    .tab-pane.active { display: block; }
    .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
    .field { display: flex; flex-direction: column; gap: 4px; margin-bottom: 10px; }
    .field label { font-size: 0.76rem; color: var(--muted); }
    .field input, .field select {
      background: var(--bg);
      border: 1px solid var(--border);
      color: var(--text);
      padding: 7px 10px;
      border-radius: 7px;
      font-size: 0.88rem;
      width: 100%;
    }
    .field input:focus, .field select:focus { outline: none; border-color: var(--cyan); }
    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(4, 1fr);
      gap: 12px;
      margin-bottom: 18px;
    }
    @media (max-width: 720px) { .kpi-grid { grid-template-columns: repeat(2, 1fr); } }
    .kpi {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 14px;
    }
    .kpi span { font-size: 0.72rem; color: var(--muted); text-transform: uppercase; }
    .kpi strong { display: block; font-size: 1.45rem; margin-top: 4px; color: #fff; }
    .kpi small { font-size: 0.75rem; color: var(--cyan); }
    .visual-grid {
      display: grid;
      grid-template-columns: 310px 1fr;
      gap: 18px;
      margin-bottom: 18px;
    }
    @media (max-width: 820px) { .visual-grid { grid-template-columns: 1fr; } }
    .svg-box {
      background: var(--bg);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 10px;
      display: flex;
      flex-direction: column;
      align-items: center;
    }
    .table-res { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
    .table-res td { padding: 7px 6px; border-bottom: 1px solid var(--border); }
    .table-res td:last-child { text-align: right; font-weight: 700; color: #fff; font-variant-numeric: tabular-nums; }
    .status-list { display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-top: 12px; }
    @media (max-width: 720px) { .status-list { grid-template-columns: 1fr; } }
    .pill {
      padding: 10px 12px;
      border-radius: 9px;
      font-size: 0.8rem;
      font-weight: 600;
      border: 1px solid;
    }
    .pill.ok { background: rgba(16, 185, 129, 0.12); border-color: var(--emerald); color: #6ee7b7; }
    .pill.warn { background: rgba(244, 63, 94, 0.14); border-color: var(--rose); color: #fda4af; }
    .error-banner {
      background: rgba(244, 63, 94, 0.18);
      border: 1px solid var(--rose);
      color: #fecdd3;
      padding: 12px;
      border-radius: 10px;
      margin-bottom: 16px;
      display: none;
      font-weight: 600;
    }
    .eq-box {
      background: var(--bg);
      border-left: 3px solid var(--amber);
      padding: 10px 12px;
      border-radius: 6px;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 0.8rem;
      margin-bottom: 8px;
      color: #3c4a48;
    }
  </style>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,300..600;1,9..144,300..500&family=JetBrains+Mono:wght@400;500;600&family=Spline+Sans:wght@400;500;600&display=swap" rel="stylesheet">
<style id="amber-editorial-theme">
  :root {
    --ink: #16201f;
    --abyss: #0e2a29;
    --teal: #1f6b6a;
    --teal-bright: #2d8b8b;
    --seafoam: #9ec9c6;
    --parchment: #f4f1e8;
    --parchment-deep: #ebe6d6;
    --cream: #fbfaf4;
    --gold: #b08440;
    --gold-soft: #c89a52;
    --rule: #d8d2bf;
    --code-bg: #10302e;
    --code-ink: #dfeae7;
    --bg: #f4f1e8;
    --card: #fbfaf4;
    --border: #d8d2bf;
    --text: #16201f;
    --muted: #5a6b69;
    --accent: #1f6b6a;
    --accent2: #b08440;
  }
  body {
    font-family: 'Spline Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif !important;
    background-color: #f4f1e8 !important;
    color: #16201f !important;
    background-image:
      radial-gradient(circle at 12% 8%, rgba(45, 139, 139, 0.06), transparent 40%),
      radial-gradient(circle at 88% 92%, rgba(176, 132, 64, 0.06), transparent 42%) !important;
    -webkit-font-smoothing: antialiased;
  }
  h1, h2, h3, .title, .brand-title {
    font-family: 'Fraunces', Georgia, serif !important;
    color: #0e2a29 !important;
    letter-spacing: -0.015em;
  }
  h1 { font-weight: 500 !important; }
  h2, h3 { font-weight: 500 !important; }
  .card, .panel, .box, section {
    background: #fbfaf4 !important;
    border: 1px solid #d8d2bf !important;
    border-radius: 12px !important;
    box-shadow: 0 12px 28px -22px rgba(14, 42, 41, 0.35) !important;
  }
  input, select {
    font-family: 'JetBrains Mono', monospace !important;
    background: #f4f1e8 !important;
    color: #0e2a29 !important;
    border: 1px solid #d8d2bf !important;
    border-radius: 8px !important;
  }
  input:focus, select:focus {
    outline: none !important;
    border-color: #1f6b6a !important;
    box-shadow: 0 0 0 3px rgba(31, 107, 106, 0.14) !important;
  }
  label, .sub, .subtitle, p {
    color: #3c4a48 !important;
  }
  code, pre, .mono, .val, .value, .metric-val, .badge, .pill, table td, table th {
    font-family: 'JetBrains Mono', monospace !important;
  }
  .badge, .pill {
    background: rgba(31, 107, 106, 0.1) !important;
    color: #1f6b6a !important;
    border: 1px solid rgba(31, 107, 106, 0.28) !important;
  }
  button, .btn {
    font-family: 'JetBrains Mono', monospace !important;
    background: #1f6b6a !important;
    color: #fbfaf4 !important;
    border: 1px solid #1f6b6a !important;
    border-radius: 8px !important;
    letter-spacing: 0.03em;
    transition: background 0.2s ease, transform 0.1s ease;
  }
  button:hover, .btn:hover {
    background: #0e2a29 !important;
  }
  button.secondary, .btn-secondary {
    background: #fbfaf4 !important;
    color: #0e2a29 !important;
    border: 1px solid #b08440 !important;
  }
  pre, .report-box, textarea {
    background: #10302e !important;
    color: #dfeae7 !important;
    border: 1px solid rgba(158, 201, 198, 0.22) !important;
    border-radius: 10px !important;
    font-family: 'JetBrains Mono', monospace !important;
  }
  hr {
    border: 0 !important;
    border-top: 1px solid #d8d2bf !important;
  }
</style>
</head>
<body>
<div id="nucoil-top-switch" style="background:rgba(244,241,232,0.92);backdrop-filter:blur(10px);border-bottom:1px solid #d8d2bf;padding:10px 18px;display:flex;justify-content:space-between;align-items:center;font-family:'Spline Sans',system-ui,sans-serif;font-size:13px;color:#0e2a29;position:sticky;top:0;z-index:40">
  <span style="font-family:'Fraunces',serif;font-size:15px;font-weight:500">NuCoil <em style="color:#b08440;font-style:italic">· Ingeniería Térmica</em></span>
  <a id="btn-abrir-plano" href="/plano" style="background:#1f6b6a;color:#fbfaf4;text-decoration:none;padding:6px 14px;border-radius:8px;font-family:'JetBrains Mono',monospace;font-size:11.5px;letter-spacing:.04em">Plano en Sección 2D →</a>
</div>
<script>
document.addEventListener('click', function(e){
  const a = e.target.closest('#btn-abrir-plano');
  if (!a) return;
  const qs = new URLSearchParams();
  document.querySelectorAll('input[id], input[name]').forEach(el => {
    if (el.value) qs.set(el.id || el.name, el.value);
  });
  const base = location.protocol === 'file:' ? 'plano.html' : '/plano';
  a.href = base + (qs.toString() ? '?' + qs.toString() : '');
});
</script>

  <header>
    <div class="brand">
      <div class="badge-cpp">C++17 ENGINE</div>
      <div>
        <h1>BioCoil Studio · Diseño de Serpentín de Enfriamiento</h1>
        <small>Fenómenos de Transporte · Diseño Térmico de Biorreactores</small>
      </div>
    </div>
    <div class="toolbar">
      <button type="button" onclick="cargarPreset('base')">Caso Base</button>
      <button type="button" onclick="cargarPreset('doran')">Verificación Doran 9.13</button>
      <button type="button" onclick="cargarPreset('alta')">Alta Densidad Celular</button>
      <button type="button" class="primary" onclick="descargarReporte()">↓ Exportar Reporte (.txt)</button>
    </div>
  </header>

  <div class="container">
    <!-- PANEL DE ENTRADAS -->
    <div class="card">
      <h2>Parámetros de Proceso <span style="font-size:0.72rem;color:var(--emerald);">● Live</span></h2>
      <div class="tabs">
        <button type="button" class="tab-btn active" onclick="abrirTab(0)">1. Calor & Rxn</button>
        <button type="button" class="tab-btn" onclick="abrirTab(1)">2. Agitación</button>
        <button type="button" class="tab-btn" onclick="abrirTab(2)">3. Serpentín</button>
      </div>

      <!-- TAB 1 -->
      <div class="tab-pane active" id="pane-0">
        <div class="field">
          <label>Método de cálculo de entalpía (ΔHᵣ)</label>
          <select id="modoEntalpia" onchange="toggleModoEntalpia(); recalcular();">
            <option value="1">Balance estequiométrico (Calores de combustión, Doran Cap. 5)</option>
            <option value="2">Ingreso directo por mol consumido (ej. base O₂)</option>
          </select>
        </div>
        <div id="grupoEstequio">
          <div class="grid-2">
            <div class="field"><label>ΔHc° Sustrato (kJ/mol)</label><input type="number" id="dhCombSustrato" value="2805" step="any"></div>
            <div class="field"><label>ΔHc° Biomasa (kJ/C-mol)</label><input type="number" id="dhCombBiomasa" value="560" step="any"></div>
          </div>
          <div class="field"><label>Rendimiento biomasa/sustrato Y_X/S (C-mol/mol)</label><input type="number" id="rendimientoYxs" value="3.0" step="any"></div>
        </div>
        <div id="grupoDirecto" style="display:none;">
          <div class="field"><label>|ΔHᵣ| directo (kJ/mol referencia)</label><input type="number" id="dhDirecto" value="460" step="any"></div>
        </div>
        <div class="grid-2">
          <div class="field"><label>Consumo molar (mol/h)</label><input type="number" id="tasaMolar" value="15.0" step="any"></div>
          <div class="field"><label>Factor calor al caldo (0–1)</label><input type="number" id="factorCalor" value="1.0" min="0" max="1" step="any"></div>
        </div>
        <hr style="border:0;border-top:1px solid var(--border);margin:12px 0;">
        <div class="grid-2">
          <div class="field"><label>Temp. del caldo (°C)</label><input type="number" id="tCaldo" value="30" step="any"></div>
          <div class="field"><label>Caudal ref. disp. (kg/s)</label><input type="number" id="flujoRefDisp" value="0.25" step="any"></div>
          <div class="field"><label>Entrada refrigerante (°C)</label><input type="number" id="tRefIn" value="10" step="any"></div>
          <div class="field"><label>Salida refrigerante (°C)</label><input type="number" id="tRefOut" value="18" step="any"></div>
        </div>
      </div>

      <!-- TAB 2 -->
      <div class="tab-pane" id="pane-1">
        <div class="grid-2">
          <div class="field"><label>Diámetro tanque Dₜ (m)</label><input type="number" id="diametroTanque" value="1.20" step="any"></div>
          <div class="field"><label>Volumen líquido Vₗ (m³)</label><input type="number" id="volumenLiquido" value="1.10" step="any"></div>
          <div class="field"><label>Diámetro impulsor Dₐ (m)</label><input type="number" id="diametroImpulsor" value="0.40" step="any"></div>
          <div class="field"><label>Velocidad agitación N (rps)</label><input type="number" id="velocidadAgitacion" value="2.5" step="any"></div>
          <div class="field"><label>Número de potencia Nₚ</label><input type="number" id="numeroPotencia" value="5.0" step="any"></div>
          <div class="field"><label>Flujo aireación Q_g (m³/s)</label><input type="number" id="flujoAire" value="0.010" step="any"></div>
        </div>
        <hr style="border:0;border-top:1px solid var(--border);margin:12px 0;">
        <div class="grid-2">
          <div class="field"><label>Densidad caldo (kg/m³)</label><input type="number" id="rhoCaldo" value="1020" step="any"></div>
          <div class="field"><label>Viscosidad caldo (Pa·s)</label><input type="number" id="muCaldo" value="0.0020" step="any"></div>
          <div class="field"><label>Conductividad caldo (W/m·K)</label><input type="number" id="kCaldo" value="0.60" step="any"></div>
          <div class="field"><label>Cp caldo (J/kg·K)</label><input type="number" id="cpCaldo" value="4100" step="any"></div>
        </div>
      </div>

      <!-- TAB 3 -->
      <div class="tab-pane" id="pane-2">
        <div class="grid-2">
          <div class="field"><label>Diámetro ext. tubo Dₒ (m)</label><input type="number" id="dExtTubo" value="0.025" step="any"></div>
          <div class="field"><label>Espesor de pared (m)</label><input type="number" id="espesorPared" value="0.002" step="any"></div>
          <div class="field"><label>Conductividad pared k (W/m·K)</label><input type="number" id="kPared" value="16.0" step="any"></div>
          <div class="field"><label>Ensuciamiento R_f (m²·K/W)</label><input type="number" id="factorEnsuciamiento" value="0.0002" step="any"></div>
          <div class="field"><label>Holgura radial pared (m)</label><input type="number" id="holguraPared" value="0.10" step="any"></div>
          <div class="field"><label>Paso entre vueltas p (m)</label><input type="number" id="pasoVueltas" value="0.075" step="any"></div>
        </div>
        <hr style="border:0;border-top:1px solid var(--border);margin:12px 0;">
        <div class="grid-2">
          <div class="field"><label>Cp refrigerante (J/kg·K)</label><input type="number" id="cpRef" value="4185" step="any"></div>
          <div class="field"><label>Densidad ref. (kg/m³)</label><input type="number" id="rhoRef" value="999" step="any"></div>
          <div class="field"><label>Viscosidad ref. (Pa·s)</label><input type="number" id="muRef" value="0.00117" step="any"></div>
          <div class="field"><label>Conductividad ref. (W/m·K)</label><input type="number" id="kRef" value="0.589" step="any"></div>
        </div>
      </div>
    </div>

    <!-- PANEL DE RESULTADOS Y ESQUEMA -->
    <div>
      <div id="errorBanner" class="error-banner"></div>

      <div class="kpi-grid">
        <div class="kpi"><span>Carga Térmica Total</span><strong id="kpiQ">—</strong><small id="kpiQSub">Reacción + Agitación</small></div>
        <div class="kpi"><span>Coeficiente Global Uₒ</span><strong id="kpiU">—</strong><small>W/(m²·K)</small></div>
        <div class="kpi"><span>Área Requerida Aₒ</span><strong id="kpiA">—</strong><small id="kpiL">Longitud: — m</small></div>
        <div class="kpi"><span>Vueltas del Serpentín</span><strong id="kpiN">—</strong><small id="kpiH">Altura: — m</small></div>
      </div>

      <div class="visual-grid">
        <div class="svg-box">
          <div style="font-size:0.78rem;color:var(--muted);margin-bottom:6px;font-weight:600;">ESQUEMA GEOMÉTRICO A ESCALA</div>
          <svg id="reactorSvg" width="270" height="310" viewBox="0 0 270 310"></svg>
          <div style="font-size:0.72rem;color:var(--muted);text-align:center;margin-top:4px;">
            <span style="color:#1f6b6a;">■</span> Caldo &nbsp;
            <span style="color:#b08440;">■</span> Serpentín sumergido &nbsp;
            <span style="color:#f43f5e;">■</span> Fuera de nivel
          </div>
        </div>

        <div class="card">
          <h2>Desglose de Ingeniería</h2>
          <table class="table-res">
            <tr><td>Entalpía de reacción calculada (|ΔHᵣ|)</td><td id="resDh">—</td></tr>
            <tr><td>Calor metabólico de fermentación (Q_rxn)</td><td id="resQrxn">—</td></tr>
            <tr><td>Reynolds de agitación (Reₐ) / Potencia aireada (P_g)</td><td id="resPg">—</td></tr>
            <tr><td>Diferencia media logarítmica (ΔT_lm)</td><td id="resDtlm">—</td></tr>
            <tr><td>Coeficiente convectivo lado caldo (hₒ)</td><td id="resHo">—</td></tr>
            <tr><td>Coeficiente convectivo interior (hᵢ) · Re tubo</td><td id="resHi">—</td></tr>
            <tr><td>Diámetro de hélice (D_hélice) / Altura líquido (H_L)</td><td id="resGeom">—</td></tr>
            <tr><td>Caudal requerido refrigerante (ṁ_req)</td><td id="resFlow">—</td></tr>
            <tr><td>Caída de presión estimada (ΔP, corr. Ito)</td><td id="resDp">—</td></tr>
          </table>

          <div class="status-list">
            <div class="pill ok" id="pillLiquid">Sumergencia: —</div>
            <div class="pill ok" id="pillImpeller">Impulsor: —</div>
            <div class="pill ok" id="pillFlow">Caudal: —</div>
          </div>
        </div>
      </div>

      <div class="card">
        <h2>Fundamentación Matemática del Modelo</h2>
        <div class="eq-box">1. Entalpía y Balance: ΔHᵣ = ΔHc°_S − Y_X/S · ΔHc°_X  |  Q_rxn = (|ΔHᵣ| · ṅ · f_calor)/3600  |  P₀ = Nₚ·ρ·N³·Dₐ⁵  |  P_g = 0.706·(P₀²·N·Dₐ³/Q_g^0.56)^0.45</div>
        <div class="eq-box">2. Convección Exterior e Interior: Nu_o = 0.87·Reₐ^0.62·Pr^(1/3) (Doran Fig. 9.17)  |  Nu_i = Nu_Gnielinski · (1 + 3.5·Dᵢ/D_hélice)</div>
        <div class="eq-box">3. Resistencias y Geometría: 1/Uₒ = 1/hₒ + R_f + [Dₒ·ln(Dₒ/Dᵢ)]/(2k) + Dₒ/(Dᵢ·hᵢ)  |  Aₒ = Q_total/(Uₒ·ΔT_lm)  |  N = L / √[(π·D_hélice)² + p²]</div>
      </div>
    </div>
  </div>

  <script>
    const campos = [
      "modoEntalpia","dhCombSustrato","dhCombBiomasa","rendimientoYxs","dhDirecto","tasaMolar","factorCalor",
      "diametroTanque","volumenLiquido","diametroImpulsor","velocidadAgitacion","numeroPotencia","flujoAire",
      "rhoCaldo","muCaldo","kCaldo","cpCaldo","tCaldo","tRefIn","tRefOut","flujoRefDisp","cpRef","rhoRef",
      "muRef","kRef","dExtTubo","espesorPared","kPared","factorEnsuciamiento","holguraPared","pasoVueltas"
    ];

    let ultimoJson = null;

    function abrirTab(idx) {
      document.querySelectorAll('.tab-btn').forEach((b, i) => b.classList.toggle('active', i === idx));
      document.querySelectorAll('.tab-pane').forEach((p, i) => p.classList.toggle('active', i === idx));
    }

    function toggleModoEntalpia() {
      const m = document.getElementById('modoEntalpia').value;
      document.getElementById('grupoEstequio').style.display = (m === '1') ? 'block' : 'none';
      document.getElementById('grupoDirecto').style.display = (m === '2') ? 'block' : 'none';
    }

    function cargarPreset(tipo) {
      const set = (id, val) => { document.getElementById(id).value = val; };
      if (tipo === 'base') {
        set('modoEntalpia', '1'); set('dhCombSustrato', 2805); set('dhCombBiomasa', 560); set('rendimientoYxs', 3.0);
        set('tasaMolar', 15.0); set('diametroTanque', 1.20); set('volumenLiquido', 1.10); set('diametroImpulsor', 0.40);
        set('velocidadAgitacion', 2.5); set('flujoAire', 0.010); set('tCaldo', 30); set('tRefIn', 10); set('tRefOut', 18);
        set('dExtTubo', 0.025); set('pasoVueltas', 0.075);
      } else if (tipo === 'doran') {
        set('modoEntalpia', '2'); set('dhDirecto', 460); set('tasaMolar', 85.0);
        set('diametroTanque', 1.50); set('volumenLiquido', 2.00); set('diametroImpulsor', 0.50);
        set('velocidadAgitacion', 3.0); set('flujoAire', 0.015); set('tCaldo', 32); set('tRefIn', 10); set('tRefOut', 20);
        set('dExtTubo', 0.030); set('pasoVueltas', 0.075);
      } else if (tipo === 'alta') {
        set('modoEntalpia', '1'); set('dhCombSustrato', 2805); set('dhCombBiomasa', 560); set('rendimientoYxs', 2.8);
        set('tasaMolar', 45.0); set('diametroTanque', 1.40); set('volumenLiquido', 1.80); set('diametroImpulsor', 0.48);
        set('velocidadAgitacion', 3.5); set('flujoAire', 0.020); set('tCaldo', 35); set('tRefIn', 8); set('tRefOut', 18);
        set('dExtTubo', 0.028); set('pasoVueltas', 0.065);
      }
      toggleModoEntalpia();
      recalcular();
    }

    function dibujarEsquema(d) {
      const svg = document.getElementById('reactorSvg');
      const cx = 135, topY = 30, botY = 275;
      const tankH_px = botY - topY;
      const maxH_m = Math.max(d.alturaTanqueEst, d.alturaSerpentin * 1.15, 0.5);
      const scaleY = tankH_px / maxH_m;
      const tankW_px = Math.min(200, Math.max(110, d.diametroTanque * 115));
      const leftX = cx - tankW_px / 2;

      const liqH_px = Math.min(tankH_px, d.alturaLiquido * scaleY);
      const liqY = botY - liqH_px;

      const coilW_px = tankW_px * (d.diametroHelice / d.diametroTanque);
      const impW_px = tankW_px * (d.diametroImpulsor / d.diametroTanque);
      const impY = botY - Math.min(liqH_px * 0.25, 35);

      let coilsSvg = '';
      const nTurns = Math.min(40, Math.max(1, Math.round(d.numVueltas)));
      const pitchPx = d.pasoVueltas * scaleY;
      const startCoilY = botY - Math.min(22, liqH_px * 0.15);

      for (let i = 0; i < nTurns; i++) {
        const y = startCoilY - i * pitchPx;
        if (y < 15) break;
        const fuera = (y < liqY);
        const color = fuera ? '#f43f5e' : '#b08440';
        coilsSvg += `<ellipse cx="${cx}" cy="${y.toFixed(1)}" rx="${(coilW_px/2).toFixed(1)}" ry="4.5" fill="none" stroke="${color}" stroke-width="3"/>`;
      }

      svg.innerHTML = `
        <rect x="${leftX}" y="${liqY.toFixed(1)}" width="${tankW_px}" height="${liqH_px.toFixed(1)}" fill="rgba(56, 189, 248, 0.18)"/>
        <line x1="${leftX}" y1="${liqY.toFixed(1)}" x2="${leftX + tankW_px}" y2="${liqY.toFixed(1)}" stroke="#1f6b6a" stroke-width="2" stroke-dasharray="4 3"/>
        <path d="M ${leftX} ${topY} L ${leftX} ${botY} Q ${cx} ${botY+16} ${leftX+tankW_px} ${botY} L ${leftX+tankW_px} ${topY}" fill="none" stroke="#6e7d7b" stroke-width="3"/>
        <line x1="${cx}" y1="${topY - 12}" x2="${cx}" y2="${impY}" stroke="#5a6b69" stroke-width="3"/>
        <rect x="${cx - impW_px/2}" y="${impY - 5}" width="${impW_px}" height="10" rx="2" fill="#3c4a48"/>
        ${coilsSvg}
        <text x="${leftX + 6}" y="${liqY - 5}" fill="#1f6b6a" font-size="10" font-family="sans-serif">Nivel caldo: ${d.alturaLiquido.toFixed(2)} m</text>
      `;
    }

    async function recalcular() {
      const params = new URLSearchParams();
      campos.forEach(id => params.append(id, document.getElementById(id).value));
      try {
        const resp = await fetch('/api/calcular?' + params.toString());
        const d = await resp.json();
        ultimoJson = d;

        const err = document.getElementById('errorBanner');
        if (!d.valido) {
          err.style.display = 'block';
          err.textContent = d.mensajeError;
          return;
        }
        err.style.display = 'none';

        document.getElementById('kpiQ').textContent = d.qTotalKw.toFixed(2) + ' kW';
        document.getElementById('kpiQSub').textContent = `Rxn: ${d.qReaccionKw.toFixed(2)} kW | Agit: ${d.potenciaAireadaKw.toFixed(2)} kW`;
        document.getElementById('kpiU').textContent = d.uGlobal.toFixed(1);
        document.getElementById('kpiA').textContent = d.areaReq.toFixed(3) + ' m²';
        document.getElementById('kpiL').textContent = `Longitud tubo: ${d.longitudTubo.toFixed(2)} m`;
        document.getElementById('kpiN').textContent = d.numVueltas.toFixed(1);
        document.getElementById('kpiH').textContent = `Altura: ${d.alturaSerpentin.toFixed(2)} / ${d.alturaLiquido.toFixed(2)} m`;

        document.getElementById('resDh').textContent = `${d.dhReaccion.toFixed(1)} kJ/mol`;
        document.getElementById('resQrxn').textContent = `${d.qReaccionKw.toFixed(2)} kW`;
        document.getElementById('resPg').textContent = `${Math.round(d.reAgitacion).toLocaleString()} / ${d.potenciaAireadaKw.toFixed(2)} kW`;
        document.getElementById('resDtlm').textContent = `${d.dtlm.toFixed(2)} K`;
        document.getElementById('resHo').textContent = `${d.ho.toFixed(1)} W/(m²·K)`;
        document.getElementById('resHi').textContent = `${d.hi.toFixed(1)} W/(m²·K) · Re=${Math.round(d.reRef).toLocaleString()}`;
        document.getElementById('resGeom').textContent = `${d.diametroHelice.toFixed(3)} m / ${d.alturaLiquido.toFixed(3)} m`;
        document.getElementById('resFlow').textContent = `${d.flujoRefReq.toFixed(3)} kg/s`;
        document.getElementById('resDp').textContent = `${d.caidaPresion.toFixed(2)} kPa`;

        const setPill = (id, ok, txtOk, txtWarn) => {
          const el = document.getElementById(id);
          el.className = 'pill ' + (ok ? 'ok' : 'warn');
          el.textContent = ok ? txtOk : txtWarn;
        };
        setPill('pillLiquid', d.cabeEnLiquido, '✓ Sumergido en caldo', '⚠ Rebasa nivel de caldo');
        setPill('pillImpeller', d.libreImpulsor, '✓ Libre de impulsor', '⚠ Muy cerca del impulsor');
        setPill('pillFlow', d.caudalSuficiente, '✓ Caudal suficiente', '⚠ Caudal insuficiente');

        dibujarEsquema(d);
      } catch (e) {
        console.error(e);
      }
    }

    function descargarReporte() {
      if (!ultimoJson || !ultimoJson.valido) return;
      const d = ultimoJson;
      const texto = [
        "====================================================================",
        "  REPORTE TECNICO DE DISENO DE SERPENTIN DE ENFRIAMIENTO (C++17)    ",
        "====================================================================",
        `Entalpia de reaccion (|dHr|)     : ${d.dhReaccion.toFixed(2)} kJ/mol`,
        `Calor metabolico (Q_rxn)         : ${d.qReaccionKw.toFixed(2)} kW`,
        `Reynolds de agitacion (Re_a)     : ${d.reAgitacion.toFixed(0)}`,
        `Potencia aireada disipada (Pg)   : ${d.potenciaAireadaKw.toFixed(2)} kW`,
        `CARGA TERMICA TOTAL (Q_total)    : ${d.qTotalKw.toFixed(2)} kW`,
        "--------------------------------------------------------------------",
        `Diferencia media log. (DT_lm)    : ${d.dtlm.toFixed(2)} K`,
        `Coeficiente lado caldo (ho)      : ${d.ho.toFixed(2)} W/(m^2*K)`,
        `Coeficiente lado tubo (hi)       : ${d.hi.toFixed(2)} W/(m^2*K)`,
        `COEFICIENTE GLOBAL (Uo)          : ${d.uGlobal.toFixed(2)} W/(m^2*K)`,
        "--------------------------------------------------------------------",
        `Area de transferencia (Ao)       : ${d.areaReq.toFixed(3)} m^2`,
        `Longitud total de tubo (L)       : ${d.longitudTubo.toFixed(3)} m`,
        `Diametro de la helice (D_h)      : ${d.diametroHelice.toFixed(3)} m`,
        `Numero de vueltas (N)            : ${d.numVueltas.toFixed(2)} vueltas`,
        `Altura serpentin / liquido       : ${d.alturaSerpentin.toFixed(3)} m / ${d.alturaLiquido.toFixed(3)} m`,
        `Caudal requerido refrigerante    : ${d.flujoRefReq.toFixed(3)} kg/s`,
        `Caida de presion estimada (dP)   : ${d.caidaPresion.toFixed(2)} kPa`,
        "===================================================================="
      ].join("\n");

      const blob = new Blob([texto], { type: "text/plain;charset=utf-8" });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = "reporte_serpentin_cpp.txt";
      a.click();
    }

    campos.forEach(id => document.getElementById(id).addEventListener('input', recalcular));
    window.addEventListener('DOMContentLoaded', recalcular);
  </script>
</body>
</html>
)HTML";

int main() {
    std::signal(SIGPIPE, SIG_IGN);

    int serverFd = -1;
    int puertoActivo = 8080;
    const char* envPort = std::getenv("PORT");
    int puertoInicial = envPort ? std::stoi(envPort) : 8080;
    int puertoFinal = envPort ? puertoInicial : 8085;

    for (int puerto = puertoInicial; puerto <= puertoFinal; ++puerto) {
        int fd = socket(AF_INET, SOCK_STREAM, 0);
        if (fd < 0) continue;

        int opt = 1;
        setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = htonl(INADDR_ANY);
        addr.sin_port = htons(puerto);

        if (bind(fd, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == 0 && listen(fd, 16) == 0) {
            serverFd = fd;
            puertoActivo = puerto;
            break;
        }
        close(fd);
    }

    if (serverFd < 0) {
        std::cerr << "Error: No se pudo iniciar el servidor en el puerto " << puertoInicial << "\n";
        return 1;
    }

    std::string url = "http://localhost:" + std::to_string(puertoActivo);
    std::cout << "====================================================================\n";
    std::cout << "  SERVIDOR NATIVO C++17 ACTIVO EN PUERTO: " << puertoActivo << "\n";
    std::cout << "====================================================================\n";

    if (!envPort) {
        std::string cmdAbrir = "open " + url + " >/dev/null 2>&1 &";
        std::system(cmdAbrir.c_str());
    }

    while (true) {
        int clientFd = accept(serverFd, nullptr, nullptr);
        if (clientFd < 0) continue;

        char buffer[8192];
        ssize_t bytesRead = recv(clientFd, buffer, sizeof(buffer) - 1, 0);
        if (bytesRead > 0) {
            buffer[bytesRead] = '\0';
            std::string req(buffer);
            std::istringstream reqLine(req);
            std::string method, target;
            reqLine >> method >> target;

            std::string body;
            std::string contentType;

            if (target.rfind("/plano", 0) == 0) {
                body = PAGINA_PLANO;
                contentType = "text/html; charset=utf-8";
            } else if (target.rfind("/api/calcular", 0) == 0) {
                ParametrosProceso p = parsearParametros(target);
                ResultadosDiseno r = calcularDiseno(p);
                body = resultadosAJson(p, r);
                contentType = "application/json; charset=utf-8";
            } else {
                body = PAGINA_HTML;
                contentType = "text/html; charset=utf-8";
            }

            std::ostringstream resp;
            resp << "HTTP/1.1 200 OK\r\n"
                 << "Content-Type: " << contentType << "\r\n"
                 << "Content-Length: " << body.size() << "\r\n"
                 << "Cache-Control: no-store\r\n"
                 << "Connection: close\r\n\r\n"
                 << body;

            std::string fullResp = resp.str();
            send(clientFd, fullResp.data(), fullResp.size(), 0);
        }
        close(clientFd);
    }

    close(serverFd);
    return 0;
}
