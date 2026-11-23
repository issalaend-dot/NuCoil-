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

const std::string PAGINA_HTML = R"HTML(<!doctype html>
<html lang="es">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BioCoil C++ | Ingeniería de Serpentines en Biorreactores</title>
  <style>
    :root {
      --bg: #0b1120;
      --panel: #131c31;
      --panel-alt: #19243e;
      --border: #263554;
      --text: #e2e8f0;
      --muted: #94a3b8;
      --cyan: #38bdf8;
      --amber: #f59e0b;
      --emerald: #10b981;
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
      background: #0f172a;
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
      background: linear-gradient(135deg, #0284c7, #2563eb);
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
    button.primary { background: #0284c7; border-color: #0284c7; color: #fff; }
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
    .tab-btn.active { background: #0284c7; border-color: #38bdf8; color: #fff; }
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
      color: #cbd5e1;
    }
  </style>
</head>
<body>
  <header>
    <div class="brand">
      <div class="badge-cpp">C++17 ENGINE</div>
      <div>
        <h1>BioCoil Studio · Diseño de Serpentín de Enfriamiento</h1>
        <small>Fenómenos de Transporte · Servidor HTTP Nativo en C++ (Sin APIs externas)</small>
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
      <h2>Parámetros de Proceso <span style="font-size:0.72rem;color:var(--emerald);">● Cálculo C++ en vivo</span></h2>
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
            <span style="color:#38bdf8;">■</span> Caldo &nbsp;
            <span style="color:#f59e0b;">■</span> Serpentín sumergido &nbsp;
            <span style="color:#f43f5e;">■</span> Fuera de nivel
          </div>
        </div>

        <div class="card">
          <h2>Desglose de Ingeniería (Resultados C++)</h2>
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
        <h2>Fundamentación Matemática Implementada en el Servidor C++</h2>
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
        const color = fuera ? '#f43f5e' : '#f59e0b';
        coilsSvg += `<ellipse cx="${cx}" cy="${y.toFixed(1)}" rx="${(coilW_px/2).toFixed(1)}" ry="4.5" fill="none" stroke="${color}" stroke-width="3"/>`;
      }

      svg.innerHTML = `
        <rect x="${leftX}" y="${liqY.toFixed(1)}" width="${tankW_px}" height="${liqH_px.toFixed(1)}" fill="rgba(56, 189, 248, 0.18)"/>
        <line x1="${leftX}" y1="${liqY.toFixed(1)}" x2="${leftX + tankW_px}" y2="${liqY.toFixed(1)}" stroke="#38bdf8" stroke-width="2" stroke-dasharray="4 3"/>
        <path d="M ${leftX} ${topY} L ${leftX} ${botY} Q ${cx} ${botY+16} ${leftX+tankW_px} ${botY} L ${leftX+tankW_px} ${topY}" fill="none" stroke="#64748b" stroke-width="3"/>
        <line x1="${cx}" y1="${topY - 12}" x2="${cx}" y2="${impY}" stroke="#94a3b8" stroke-width="3"/>
        <rect x="${cx - impW_px/2}" y="${impY - 5}" width="${impW_px}" height="10" rx="2" fill="#cbd5e1"/>
        ${coilsSvg}
        <text x="${leftX + 6}" y="${liqY - 5}" fill="#38bdf8" font-size="10" font-family="sans-serif">Nivel caldo: ${d.alturaLiquido.toFixed(2)} m</text>
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
</html>)HTML";

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

            if (target.rfind("/api/calcular", 0) == 0) {
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
