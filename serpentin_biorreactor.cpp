#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <limits>

constexpr double PI = 3.14159265358979323846;

struct ParametrosProceso {
    // 1. Estequiometria y generacion de calor (Doran Cap. 5 y Cap. 8)
    int modoEntalpia = 1;            // 1: Calculo estequiometrico (combustiones), 2: Directo
    double dhCombSustrato = 2805.0;  // kJ/mol sustrato (ej. glucosa = 2805)
    double dhCombBiomasa = 560.0;    // kJ/C-mol biomasa (ej. CH1.8O0.5N0.2 = 560)
    double rendimientoYxs = 3.0;     // C-mol biomasa / mol sustrato
    double dhDirecto = 460.0;        // kJ/mol referencia (si modoEntalpia == 2)
    double tasaMolar = 15.0;         // mol/h de sustrato (o especie de referencia)
    double factorCalor = 1.0;        // Fraccion de calor metabolico liberado al caldo

    // 2. Agitacion y propiedades del caldo (Transferencia de momento - Doran Cap. 8)
    double diametroTanque = 1.20;    // m (Dt)
    double volumenLiquido = 1.10;    // m^3 (VL)
    double diametroImpulsor = 0.40;  // m (Da, tipicamente Dt/3)
    double velocidadAgitacion = 2.5; // rps (150 rpm)
    double numeroPotencia = 5.0;     // Np (5.0 para turbina Rushton turbulenta)
    double flujoAire = 0.010;        // m^3/s (Qg, 0 si no hay aireacion)
    double rhoCaldo = 1020.0;        // kg/m^3
    double muCaldo = 0.0020;         // Pa*s (2.0 mPa*s)
    double kCaldo = 0.60;            // W/(m*K)
    double cpCaldo = 4100.0;         // J/(kg*K)

    // 3. Temperaturas y propiedades del refrigerante
    double tCaldo = 30.0;            // °C
    double tRefIn = 10.0;            // °C
    double tRefOut = 18.0;           // °C
    double flujoRefDisp = 0.25;      // kg/s disponible
    double cpRef = 4185.0;           // J/(kg*K)
    double rhoRef = 999.0;           // kg/m^3
    double muRef = 0.00117;          // Pa*s (1.17 mPa*s)
    double kRef = 0.589;             // W/(m*K)

    // 4. Geometria y material del serpentin
    double dExtTubo = 0.025;         // m (25 mm)
    double espesorPared = 0.002;     // m (2 mm)
    double kPared = 16.0;            // W/(m*K) Acero inoxidable 316
    double factorEnsuciamiento = 0.0002; // m^2*K/W
    double holguraPared = 0.10;      // m (100 mm)
    double pasoVueltas = 0.075;      // m (75 mm)
};

struct ResultadosDiseno {
    bool valido = true;
    std::string mensajeError;

    // Balance energetico
    double dhReaccion = 0.0;     // kJ/mol
    double qReaccion = 0.0;      // W
    double reAgitacion = 0.0;    // adimensional
    double potenciaNoAir = 0.0;  // W (P0)
    double potenciaAireada = 0.0;// W (Pg)
    double qTotal = 0.0;         // W

    // Transferencia de calor
    double dt1 = 0.0;            // K
    double dt2 = 0.0;            // K
    double dtlm = 0.0;           // K
    double prCaldo = 0.0;
    double nuCaldo = 0.0;
    double ho = 0.0;             // W/(m^2*K)
    double flujoRefReq = 0.0;    // kg/s
    double velRef = 0.0;         // m/s
    double reRef = 0.0;
    double prRef = 0.0;
    double hi = 0.0;             // W/(m^2*K)
    double uGlobal = 0.0;        // W/(m^2*K)

    // Geometria e hidraulica
    double areaReq = 0.0;        // m^2
    double longitudTubo = 0.0;   // m
    double diametroHelice = 0.0; // m
    double numVueltas = 0.0;
    double alturaSerpentin = 0.0;// m
    double alturaLiquido = 0.0;  // m
    double caidaPresion = 0.0;   // kPa
    bool cabeEnLiquido = true;
    bool libreImpulsor = true;
    bool caudalSuficiente = true;
};

ResultadosDiseno calcularDiseno(const ParametrosProceso& p) {
    ResultadosDiseno r;

    // Validaciones fisicas basicas
    if (p.tCaldo <= p.tRefIn || p.tCaldo <= p.tRefOut || p.tRefOut <= p.tRefIn) {
        r.valido = false;
        r.mensajeError = "Error termodinamico: Se requiere T_caldo > T_salida > T_entrada.";
        return r;
    }
    double dIntTubo = p.dExtTubo - 2.0 * p.espesorPared;
    if (dIntTubo <= 0.0) {
        r.valido = false;
        r.mensajeError = "Error geometrico: El espesor de pared supera el radio del tubo.";
        return r;
    }

    // --- 1. BALANCE DE ENERGIA Y ENTALPIA DE REACCION ---
    if (p.modoEntalpia == 1) {
        // Balance por calores de combustion (Doran Cap. 5): dHr = dHc_S - Y_X/S * dHc_X
        r.dhReaccion = p.dhCombSustrato - (p.rendimientoYxs * p.dhCombBiomasa);
    } else {
        r.dhReaccion = p.dhDirecto;
    }
    if (r.dhReaccion < 0.0) r.dhReaccion = 0.0;

    // Q_rxn en Watts: (kJ/mol * 1000 J/kJ * mol/h) / 3600 s/h
    r.qReaccion = (r.dhReaccion * 1000.0 * p.tasaMolar * p.factorCalor) / 3600.0;

    // Potencia de agitacion (Transferencia de momento - Doran Cap. 8)
    r.reAgitacion = (p.rhoCaldo * p.velocidadAgitacion * std::pow(p.diametroImpulsor, 2)) / p.muCaldo;
    r.potenciaNoAir = p.numeroPotencia * p.rhoCaldo * std::pow(p.velocidadAgitacion, 3) * std::pow(p.diametroImpulsor, 5);

    if (p.flujoAire > 0.0) {
        // Correlacion de Michel-Miller para potencia aireada Pg
        double term = (std::pow(r.potenciaNoAir, 2) * p.velocidadAgitacion * std::pow(p.diametroImpulsor, 3))
                      / std::pow(p.flujoAire, 0.56);
        r.potenciaAireada = 0.706 * std::pow(term, 0.45);
        if (r.potenciaAireada > r.potenciaNoAir) {
            r.potenciaAireada = r.potenciaNoAir;
        }
    } else {
        r.potenciaAireada = r.potenciaNoAir;
    }

    r.qTotal = r.qReaccion + r.potenciaAireada;

    // --- 2. FUERZA IMPULSORA (DTLM) ---
    r.dt1 = p.tCaldo - p.tRefIn;
    r.dt2 = p.tCaldo - p.tRefOut;
    if (std::abs(r.dt1 - r.dt2) < 1e-6) {
        r.dtlm = 0.5 * (r.dt1 + r.dt2);
    } else {
        r.dtlm = (r.dt1 - r.dt2) / std::log(r.dt1 / r.dt2);
    }

    // --- 3. COEFICIENTES DE TRANSFERENCIA DE CALOR ---
    // Lado caldo (exterior del serpentin): Correlacion de Doran Fig. 9.17 / Problema 9.13
    r.prCaldo = (p.cpCaldo * p.muCaldo) / p.kCaldo;
    r.nuCaldo = 0.87 * std::pow(r.reAgitacion, 0.62) * std::pow(r.prCaldo, 1.0 / 3.0);
    r.ho = (r.nuCaldo * p.kCaldo) / p.diametroTanque;

    // Lado refrigerante (interior del tubo helicoidal)
    r.flujoRefReq = r.qTotal / (p.cpRef * (p.tRefOut - p.tRefIn));
    r.caudalSuficiente = (p.flujoRefDisp >= r.flujoRefReq);

    double areaTransversal = PI * std::pow(dIntTubo, 2) / 4.0;
    r.velRef = r.flujoRefReq / (p.rhoRef * areaTransversal);
    r.reRef = (p.rhoRef * r.velRef * dIntTubo) / p.muRef;
    r.prRef = (p.cpRef * p.muRef) / p.kRef;

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
        // Gnielinski + correccion por curvatura helicoidal
        fDarcy = std::pow(0.79 * std::log(r.reRef) - 1.64, -2.0);
        double num = (fDarcy / 8.0) * (r.reRef - 1000.0) * r.prRef;
        double den = 1.0 + 12.7 * std::sqrt(fDarcy / 8.0) * (std::pow(r.prRef, 2.0 / 3.0) - 1.0);
        double nuRecto = num / den;
        double factorHelice = 1.0 + 3.5 * (dIntTubo / r.diametroHelice);
        nuInt = nuRecto * factorHelice;
    }
    r.hi = (nuInt * p.kRef) / dIntTubo;

    // Coeficiente global Uo referido al area exterior
    double resExt = 1.0 / r.ho;
    double resPared = (p.dExtTubo * std::log(p.dExtTubo / dIntTubo)) / (2.0 * p.kPared);
    double resInt = p.dExtTubo / (dIntTubo * r.hi);
    double rTotal = resExt + p.factorEnsuciamiento + resPared + resInt;
    r.uGlobal = 1.0 / rTotal;

    // --- 4. DIMENSIONAMIENTO GEOMETRICO Y VALIDACION ---
    r.areaReq = r.qTotal / (r.uGlobal * r.dtlm);
    r.longitudTubo = r.areaReq / (PI * p.dExtTubo);

    double longPorVuelta = std::sqrt(std::pow(PI * r.diametroHelice, 2) + std::pow(p.pasoVueltas, 2));
    r.numVueltas = r.longitudTubo / longPorVuelta;
    r.alturaSerpentin = r.numVueltas * p.pasoVueltas;
    r.alturaLiquido = (4.0 * p.volumenLiquido) / (PI * std::pow(p.diametroTanque, 2));

    // Caida de presion con factor de curvatura helicoidal de Ito
    double factorIto = 1.0 + 0.0823 * std::pow(dIntTubo / r.diametroHelice, 0.53) * std::pow(r.reRef, 0.25);
    double dpPa = (fDarcy * factorIto) * (r.longitudTubo / dIntTubo) * (0.5 * p.rhoRef * std::pow(r.velRef, 2));
    r.caidaPresion = dpPa / 1000.0;

    r.cabeEnLiquido = (r.alturaSerpentin < r.alturaLiquido);
    r.libreImpulsor = (r.diametroHelice > p.diametroImpulsor * 1.15);

    return r;
}

void imprimirResultados(const ParametrosProceso& p, const ResultadosDiseno& r, std::ostream& out) {
    out << std::fixed << std::setprecision(2);
    out << "\n====================================================================\n";
    out << "   DISENO TERMICO E HIDRAULICO DE SERPENTIN PARA BIORREACTOR (C++)  \n";
    out << "====================================================================\n";

    if (!r.valido) {
        out << "[ERROR DE DISENO]: " << r.mensajeError << "\n";
        return;
    }

    out << "\n--- 1. BALANCE DE ENERGIA Y ESTEQUIOMETRIA ---\n";
    out << "  Entalpia de reaccion (|dHr|)     : " << r.dhReaccion << " kJ/mol\n";
    out << "  Calor metabolico (Q_rxn)         : " << r.qReaccion / 1000.0 << " kW (" << r.qReaccion << " W)\n";
    out << "  Reynolds de agitacion (Re_a)     : " << r.reAgitacion << "\n";
    out << "  Potencia sin airear (P0)         : " << r.potenciaNoAir / 1000.0 << " kW\n";
    out << "  Potencia aireada disipada (Pg)   : " << r.potenciaAireada / 1000.0 << " kW\n";
    out << "  CARGA TERMICA TOTAL (Q_total)    : " << r.qTotal / 1000.0 << " kW\n";

    out << "\n--- 2. TRANSFERENCIA DE CALOR ---\n";
    out << "  Diferencia media log. (DT_lm)    : " << r.dtlm << " K\n";
    out << "  Coeficiente lado caldo (ho)      : " << r.ho << " W/(m^2*K)\n";
    out << "  Reynolds lado refrigerante (Re)  : " << r.reRef << "\n";
    out << "  Coeficiente lado tubo (hi)       : " << r.hi << " W/(m^2*K)\n";
    out << "  COEFICIENTE GLOBAL (Uo)          : " << r.uGlobal << " W/(m^2*K)\n";

    out << "\n--- 3. GEOMETRIA DEL SERPENTIN Y VALIDACION ---\n";
    out << std::setprecision(3);
    out << "  Area de transferencia (Ao)       : " << r.areaReq << " m^2\n";
    out << "  Longitud total de tubo (L)       : " << r.longitudTubo << " m\n";
    out << "  Diametro de la helice (D_h)      : " << r.diametroHelice << " m\n";
    out << std::setprecision(2);
    out << "  Numero estimado de vueltas (N)   : " << r.numVueltas << " vueltas\n";
    out << std::setprecision(3);
    out << "  Altura ocupada por serpentin (H) : " << r.alturaSerpentin << " m (Altura liquido: " << r.alturaLiquido << " m)\n";
    out << "  Caudal requerido refrigerante    : " << r.flujoRefReq << " kg/s (Disponible: " << p.flujoRefDisp << " kg/s)\n";
    out << std::setprecision(2);
    out << "  Caida de presion estimada (dP)   : " << r.caidaPresion << " kPa\n";

    out << "\n--- 4. DIAGNOSTICO DE FACTIBILIDAD ---\n";
    out << "  [Sumergencia en caldo]  : " << (r.cabeEnLiquido ? "OK (H_serpentin < H_liquido)" : "ALERTA: El serpentin rebasa el nivel del caldo") << "\n";
    out << "  [Holgura con impulsor]  : " << (r.libreImpulsor ? "OK (Sin interferencia con impulsor)" : "ALERTA: Diametro de helice muy cercano al impulsor") << "\n";
    out << "  [Caudal refrigerante]   : " << (r.caudalSuficiente ? "OK (Caudal disponible suficiente)" : "ALERTA: Aumentar flujo o reducir T_entrada") << "\n";
    out << "====================================================================\n";
}

void exportarReporte(const ParametrosProceso& p, const ResultadosDiseno& r) {
    std::ofstream archivo("reporte_serpentin.txt");
    if (!archivo) {
        std::cout << "No se pudo crear el archivo reporte_serpentin.txt\n";
        return;
    }
    imprimirResultados(p, r, archivo);
    archivo.close();
    std::cout << "\n>> Reporte guardado exitosamente en 'reporte_serpentin.txt'\n";
}

void editarParametros(ParametrosProceso& p) {
    int sub = -1;
    while (sub != 0) {
        std::cout << "\n--- MODIFICAR CONDICIONES DE OPERACION ---\n";
        std::cout << "1. Modo de entalpia (1=Estequiometrico, 2=Directo) [Actual: " << p.modoEntalpia << "]\n";
        std::cout << "2. Tasa molar de consumo (mol/h)                   [Actual: " << p.tasaMolar << "]\n";
        std::cout << "3. Velocidad de agitacion (rps)                    [Actual: " << p.velocidadAgitacion << "]\n";
        std::cout << "4. Flujo de aireacion Qg (m^3/s)                   [Actual: " << p.flujoAire << "]\n";
        std::cout << "5. Temperatura del caldo (°C)                      [Actual: " << p.tCaldo << "]\n";
        std::cout << "6. Temperaturas refrigerante Tin / Tout (°C)       [Actual: " << p.tRefIn << " / " << p.tRefOut << "]\n";
        std::cout << "7. Diametro del tanque (m) y Volumen liquido (m^3) [Actual: " << p.diametroTanque << " / " << p.volumenLiquido << "]\n";
        std::cout << "8. Diametro exterior tubo (m) y Paso vueltas (m)   [Actual: " << p.dExtTubo << " / " << p.pasoVueltas << "]\n";
        std::cout << "0. Volver y recalcular\n";
        std::cout << "Elige una opcion: ";
        if (!(std::cin >> sub)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        switch (sub) {
            case 1:
                std::cout << "Modo (1 = por calores de combustion, 2 = directo): ";
                std::cin >> p.modoEntalpia;
                if (p.modoEntalpia == 1) {
                    std::cout << "dHc sustrato (kJ/mol), dHc biomasa (kJ/C-mol), Y_X/S (C-mol/mol): ";
                    std::cin >> p.dhCombSustrato >> p.dhCombBiomasa >> p.rendimientoYxs;
                } else {
                    std::cout << "|dHr| directo (kJ/mol): ";
                    std::cin >> p.dhDirecto;
                }
                break;
            case 2:
                std::cout << "Nueva tasa molar (mol/h): ";
                std::cin >> p.tasaMolar;
                break;
            case 3:
                std::cout << "Nueva velocidad de agitacion (rps): ";
                std::cin >> p.velocidadAgitacion;
                break;
            case 4:
                std::cout << "Nuevo flujo de aireacion Qg (m^3/s): ";
                std::cin >> p.flujoAire;
                break;
            case 5:
                std::cout << "Nueva temperatura del caldo (°C): ";
                std::cin >> p.tCaldo;
                break;
            case 6:
                std::cout << "Nueva T_entrada y T_salida del refrigerante (°C): ";
                std::cin >> p.tRefIn >> p.tRefOut;
                break;
            case 7:
                std::cout << "Nuevo diametro de tanque (m) y volumen liquido (m^3): ";
                std::cin >> p.diametroTanque >> p.volumenLiquido;
                break;
            case 8:
                std::cout << "Nuevo diametro exterior tubo (m) y paso entre vueltas (m): ";
                std::cin >> p.dExtTubo >> p.pasoVueltas;
                break;
            default:
                break;
        }
    }
}

int main() {
    ParametrosProceso params;
    int opcion = -1;

    while (opcion != 0) {
        ResultadosDiseno res = calcularDiseno(params);
        imprimirResultados(params, res, std::cout);

        std::cout << "\nMENU PRINCIPAL:\n";
        std::cout << "  [1] Modificar condiciones de operacion\n";
        std::cout << "  [2] Cargar caso de verificacion (Fermentacion aerobica tipo Doran 9.13)\n";
        std::cout << "  [3] Exportar reporte tecnico (.txt)\n";
        std::cout << "  [0] Salir\n";
        std::cout << "Selecciona una opcion: ";

        if (!(std::cin >> opcion)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        if (opcion == 1) {
            editarParametros(params);
        } else if (opcion == 2) {
            params = ParametrosProceso(); // Restaura caso base estandarizado
            params.diametroTanque = 1.50;
            params.volumenLiquido = 2.00;
            params.diametroImpulsor = 0.50;
            params.velocidadAgitacion = 3.0;
            params.tasaMolar = 28.0;
            std::cout << "\n>> Caso de verificacion cargado.\n";
        } else if (opcion == 3) {
            exportarReporte(params, res);
        }
    }

    return 0;
}