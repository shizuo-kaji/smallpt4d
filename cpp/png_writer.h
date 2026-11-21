#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace simple_png {

inline void append_be_u32(std::vector<uint8_t>& out, uint32_t value) {
    out.push_back(static_cast<uint8_t>((value >> 24) & 0xFFu));
    out.push_back(static_cast<uint8_t>((value >> 16) & 0xFFu));
    out.push_back(static_cast<uint8_t>((value >> 8) & 0xFFu));
    out.push_back(static_cast<uint8_t>(value & 0xFFu));
}

inline uint32_t crc32(const uint8_t* data, std::size_t size) {
    uint32_t crc = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) {
            if (crc & 1u) {
                crc = (crc >> 1u) ^ 0xEDB88320u;
            } else {
                crc >>= 1u;
            }
        }
    }
    return crc ^ 0xFFFFFFFFu;
}

inline uint32_t adler32(const uint8_t* data, std::size_t size) {
    const uint32_t mod = 65521u;
    uint32_t a = 1u;
    uint32_t b = 0u;
    for (std::size_t i = 0; i < size; ++i) {
        a += data[i];
        if (a >= mod) {
            a -= mod;
        }
        b += a;
        if (b >= mod) {
            b -= mod;
        }
    }
    return (b << 16u) | a;
}

inline void append_chunk(
    std::vector<uint8_t>& out,
    const char type[4],
    const std::vector<uint8_t>& payload
) {
    append_be_u32(out, static_cast<uint32_t>(payload.size()));
    const std::size_t type_offset = out.size();
    out.push_back(static_cast<uint8_t>(type[0]));
    out.push_back(static_cast<uint8_t>(type[1]));
    out.push_back(static_cast<uint8_t>(type[2]));
    out.push_back(static_cast<uint8_t>(type[3]));
    out.insert(out.end(), payload.begin(), payload.end());
    const uint32_t crc = crc32(out.data() + type_offset, 4u + payload.size());
    append_be_u32(out, crc);
}

inline bool write_png_rgb8(
    const std::string& filename,
    int width,
    int height,
    const uint8_t* pixels,
    int stride
) {
    if (width <= 0 || height <= 0 || pixels == 0 || stride < width * 3) {
        return false;
    }

    std::vector<uint8_t> raw_data;
    raw_data.reserve(static_cast<std::size_t>(height) * (1u + static_cast<std::size_t>(width) * 3u));
    for (int y = 0; y < height; ++y) {
        raw_data.push_back(0u); // filter type 0 (None)
        const uint8_t* row = pixels + static_cast<std::size_t>(y) * static_cast<std::size_t>(stride);
        raw_data.insert(raw_data.end(), row, row + static_cast<std::size_t>(width) * 3u);
    }

    std::vector<uint8_t> zlib_stream;
    zlib_stream.reserve(raw_data.size() + raw_data.size() / 65535u * 5u + 6u);
    zlib_stream.push_back(0x78u); // zlib header
    zlib_stream.push_back(0x01u); // no preset dictionary, low compression

    std::size_t cursor = 0;
    while (cursor < raw_data.size()) {
        const std::size_t remaining = raw_data.size() - cursor;
        const uint16_t block_size = static_cast<uint16_t>(std::min<std::size_t>(remaining, 65535u));
        const bool final_block = (cursor + block_size == raw_data.size());

        zlib_stream.push_back(final_block ? 0x01u : 0x00u); // BFINAL + uncompressed block type
        zlib_stream.push_back(static_cast<uint8_t>(block_size & 0xFFu));
        zlib_stream.push_back(static_cast<uint8_t>((block_size >> 8) & 0xFFu));
        const uint16_t nlen = static_cast<uint16_t>(~block_size);
        zlib_stream.push_back(static_cast<uint8_t>(nlen & 0xFFu));
        zlib_stream.push_back(static_cast<uint8_t>((nlen >> 8) & 0xFFu));
        zlib_stream.insert(
            zlib_stream.end(),
            raw_data.begin() + static_cast<std::ptrdiff_t>(cursor),
            raw_data.begin() + static_cast<std::ptrdiff_t>(cursor + block_size)
        );
        cursor += block_size;
    }

    const uint32_t check = adler32(raw_data.data(), raw_data.size());
    append_be_u32(zlib_stream, check);

    std::vector<uint8_t> png;
    png.reserve(8u + 25u + 12u + zlib_stream.size() + 12u);
    const uint8_t signature[] = {137u, 80u, 78u, 71u, 13u, 10u, 26u, 10u};
    png.insert(png.end(), signature, signature + 8u);

    std::vector<uint8_t> ihdr;
    ihdr.reserve(13u);
    append_be_u32(ihdr, static_cast<uint32_t>(width));
    append_be_u32(ihdr, static_cast<uint32_t>(height));
    ihdr.push_back(8u); // bit depth
    ihdr.push_back(2u); // color type: RGB
    ihdr.push_back(0u); // compression method
    ihdr.push_back(0u); // filter method
    ihdr.push_back(0u); // interlace method

    append_chunk(png, "IHDR", ihdr);
    append_chunk(png, "IDAT", zlib_stream);
    append_chunk(png, "IEND", std::vector<uint8_t>());

    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file) {
        return false;
    }

    file.write(reinterpret_cast<const char*>(png.data()), static_cast<std::streamsize>(png.size()));
    return file.good();
}

} // namespace simple_png
