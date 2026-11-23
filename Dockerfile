FROM gcc:13-bookworm
WORKDIR /app
COPY serpentin_web.cpp .
RUN g++ -std=c++17 -O2 serpentin_web.cpp -o serpentin_web
EXPOSE 8080
CMD ["./serpentin_web"]
